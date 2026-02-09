#!/usr/bin/env bash
# ============================================================================
# 1-orca-init-opt.sh — Gas-phase geometry optimisation (Stage 1)
# ============================================================================
#
# OVERVIEW
#   For each supplied molecule this script:
#     1. Reads charge and multiplicity from the XYZ comment line,
#     2. Writes an ORCA 6 geometry-optimisation input file,
#     3. Either runs ORCA directly (--local) or submits a SLURM job.
#
#   Execution mode is chosen automatically: if sbatch is available and
#   --local is not set, SLURM mode is used.  Otherwise, ORCA is run
#   directly in the current shell.
#
# Usage:
#   1-orca-init-opt.sh TAG
#   1-orca-init-opt.sh path/to/foo.xyz --method PBE0 --grid 2
#   1-orca-init-opt.sh --list molecules.txt --dry-run
#   1-orca-init-opt.sh --local --cpus 4 aspirin
#
# Flags:
#   -m | --method METHOD           DFT functional               [B3LYP]
#   -b | --basis  BASIS            Basis set                    [def2-SVP def2/J]
#        --disp {auto,none,D3BJ}   Dispersion correction        [auto]
#        --max-iter N              SCF iteration limit           [300]
#   -c | --cpus N                  CPU cores (%pal + SLURM)     [4]
#   -g | --grid N                  DEFGRID level (1–3)          [3]
#        --mem-per-cpu MB          Memory per core in MB         [2048]
#        --orca-bin PATH           Path to ORCA binary           [auto]
#        --list FILE               File of TAGs or XYZ paths
#        --local                   Run ORCA directly (no SLURM)
#        --dry-run                 Write inputs but do not run
#   -h | --help                    Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --partition NAME          SLURM partition               [circe]
#        --time HH:MM:SS           Wall-clock limit              [06:00:00]
#
# XYZ charge/multiplicity convention:
#   Line 2 of each .xyz file may specify charge and multiplicity in
#   either of two formats:
#     charge=0 mult=1       (key=value tokens, any order)
#     0 1                   (bare integers: charge multiplicity)
#   If neither is found, the defaults charge=0 mult=1 are used.
#
# Directory layout produced:
#   <TAG>/
#   └── <TAG>_orca_opt/
#       ├── <TAG>.inp          ORCA input
#       ├── <TAG>.log          ORCA output  (after execution)
#       ├── <TAG>.xyz          optimised geometry (ORCA writes this)
#       └── <TAG>.slurm        SLURM script (HPC mode only)
#
# Examples:
#   # Local workstation, 4 cores, default method
#   1-orca-init-opt.sh --local --cpus 4 aspirin
#
#   # HPC cluster, PBE0 functional, submit to SLURM
#   1-orca-init-opt.sh --method PBE0 --partition gpu aspirin
#
#   # Dry run: inspect inputs without running anything
#   1-orca-init-opt.sh --dry-run --local --list molecules.txt
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_METHOD="B3LYP"
DEFAULT_BASIS="def2-SVP def2/J"
DEFAULT_DISP="auto"
DEFAULT_MAX_ITER=300
DEFAULT_CPUS=4
DEFAULT_GRID=3
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="circe"
DEFAULT_WALL="06:00:00"

XYZ_DIR="pre_xyz"           # where starting geometries live

# ============================================================================
# ORCA BINARY RESOLUTION
# ============================================================================
# Priority: --orca-bin flag  >  ORCA_BIN env var  >  `which orca`
# The variable orca_bin is set during CLI parsing and finalised in main().
resolve_orca_bin() {
    # $1 = value from --orca-bin flag (may be empty)
    local flag_val=$1

    if [[ -n $flag_val ]]; then
        printf '%s' "$flag_val"
        return
    fi

    if [[ -n ${ORCA_BIN:-} ]]; then
        printf '%s' "$ORCA_BIN"
        return
    fi

    # fall back to PATH lookup
    local found
    found=$(command -v orca 2>/dev/null || true)
    if [[ -n $found ]]; then
        printf '%s' "$found"
        return
    fi

    # nothing found — caller will decide whether this is fatal
    printf ''
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
die()  { printf 'Error: %s\n' "$*" >&2; exit 1; }
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
warn() { printf '[%s] Warning: %s\n' "$(date '+%F %T')" "$*" >&2; }

require_file() { [[ -f $1 ]] || die "File '$1' not found"; }

show_help() {
    sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# \{0,1\}//'
    exit 0
}

# ============================================================================
# EXECUTION-MODE DETECTION
# ============================================================================
# Returns "local" or "slurm".  --local flag always wins; otherwise we
# check whether sbatch is available.
detect_mode() {
    local force_local=$1
    if $force_local; then
        printf 'local'
    elif command -v sbatch >/dev/null 2>&1; then
        printf 'slurm'
    else
        printf 'local'
    fi
}

# ============================================================================
# CHARGE / MULTIPLICITY PARSER
# ============================================================================
# Reads line 2 of an XYZ file and sets the global variables `charge`
# and `mult`.  Supports two formats:
#   charge=0 mult=1     (key=value, any surrounding text)
#   0 1                 (bare integers)
read_xyz_header() {
    local xyz_file=$1
    local header_line
    header_line=$(sed -n '2p' "$xyz_file")

    # try key=value format first
    if [[ $header_line =~ charge[[:space:]]*=[[:space:]]*([+-]?[0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
    else
        charge=0
    fi

    if [[ $header_line =~ mult[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
        mult=${BASH_REMATCH[1]}
    # fall back to bare "int int" format
    elif [[ $header_line =~ ^[[:space:]]*([+-]?[0-9]+)[[:space:]]+([0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
        mult=${BASH_REMATCH[2]}
    else
        mult=1
    fi
}

# ============================================================================
# DISPERSION LOGIC
# ============================================================================
# Returns the dispersion keyword to append to the ORCA bang line, or
# an empty string if none is needed.  In "auto" mode the script avoids
# doubling up on functionals that already include dispersion.
disp_keyword() {
    local method_upper=${method^^}
    case $disp_mode in
        none|NONE)
            printf '' ;;
        D3BJ|d3bj)
            printf ' D3BJ' ;;
        D4|d4)
            printf ' D4' ;;
        auto|AUTO)
            # skip if the method string already contains a dispersion suffix
            if [[ $method =~ (-D[0-9]?|-D3BJ|-D3ZERO|-D4)($|[[:space:]]) ]]; then
                printf ''; return
            fi
            # skip for functionals with built-in dispersion
            case $method_upper in
                WB97X-D|WB97X-D3|WB97X-D4|WB97X-V|WB97XD|WB97M-V|B97-D|B97-D3)
                    printf ''; return ;;
            esac
            # default: add D3BJ
            printf ' D3BJ'
            ;;
        *)
            die "--disp must be auto, none, D3BJ, or D4 (got '$disp_mode')" ;;
    esac
}

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    method=$DEFAULT_METHOD
    basis=$DEFAULT_BASIS
    disp_mode=$DEFAULT_DISP
    max_iter=$DEFAULT_MAX_ITER
    cpus=$DEFAULT_CPUS
    grid=$DEFAULT_GRID
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=$DEFAULT_PARTITION
    wall=$DEFAULT_WALL
    dry_run=false
    force_local=false
    list_file=""
    xyz_arg=""
    orca_bin_flag=""        # raw value from --orca-bin

    local opts
    opts=$(getopt -o hb:m:c:g: \
        --long help,basis:,method:,disp:,max-iter:,cpus:,grid:,\
mem-per-cpu:,partition:,time:,list:,dry-run,local,orca-bin: -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)      method=$2;         shift 2 ;;
            -b|--basis)       basis=$2;          shift 2 ;;
            --disp)           disp_mode=$2;      shift 2 ;;
            --max-iter)       max_iter=$2;       shift 2 ;;
            -c|--cpus)        cpus=$2;           shift 2 ;;
            -g|--grid)        grid=$2;           shift 2 ;;
            --mem-per-cpu)    mem_mb=$2;         shift 2 ;;
            --partition)      partition=$2;       shift 2 ;;
            --time)           wall=$2;           shift 2 ;;
            --list)           list_file=$2;      shift 2 ;;
            --orca-bin)       orca_bin_flag=$2;  shift 2 ;;
            --local)          force_local=true;  shift ;;
            --dry-run)        dry_run=true;      shift ;;
            -h|--help)        show_help ;;
            --)               shift; break ;;
            *)                die "Unknown option '$1'" ;;
        esac
    done

    # positional: either --list was given, or exactly one TAG/XYZ
    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional arguments not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide exactly one TAG or XYZ path (got $#)"
        if [[ $1 == *.xyz ]]; then
            xyz_arg=$1
        else
            xyz_arg="${XYZ_DIR}/$1.xyz"
        fi
        require_file "$xyz_arg"
    fi

    # clamp grid to 1–3
    (( grid < 1 )) && grid=1 || true
    (( grid > 3 )) && grid=3 || true
}

# ============================================================================
# ORCA INPUT WRITER
# ============================================================================
write_orca_input() {
    local inp_file=$1 xyz_file=$2
    local disp
    disp=$(disp_keyword)

    cat >"$inp_file" <<EOF
! ${method} ${basis} Opt TightSCF RIJCOSX DEFGRID${grid}${disp}

%scf
    MaxIter ${max_iter}
end

%pal nprocs ${cpus} end
%maxcore ${mem_mb}

* xyz ${charge} ${mult}
EOF
    tail -n +3 "$xyz_file" >>"$inp_file"
    echo "*" >>"$inp_file"
}

# ============================================================================
# SLURM SCRIPT WRITER  (HPC mode only)
# ============================================================================
write_slurm() {
    local slurm_file=$1 tag=$2

    # resolve the ORCA root directory for the SLURM environment block
    local orca_dir
    orca_dir=$(dirname "$orca_bin")

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=opt_${tag}
#SBATCH --partition=${partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${cpus}
#SBATCH --mem-per-cpu=${mem_mb}
#SBATCH --time=${wall}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ---- ORCA + OpenMPI environment ----
export PATH="${orca_dir}:\$PATH"
export LD_LIBRARY_PATH="${orca_dir}:\$LD_LIBRARY_PATH"

"${orca_bin}" "${tag}.inp" > "${tag}.log"
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# PROCESS ONE MOLECULE
# ============================================================================
process_molecule() {
    local xyz_file=$1 tag=$2

    local workdir="${tag}/${tag}_orca_opt"
    mkdir -p "$workdir"

    local inp_file="${workdir}/${tag}.inp"
    local log_file="${workdir}/${tag}.log"

    # parse charge / multiplicity from the XYZ comment line
    read_xyz_header "$xyz_file"

    # write the ORCA input
    write_orca_input "$inp_file" "$xyz_file"

    # ---- dry run: show what would happen and return --------------------
    if $dry_run; then
        log "[${tag}] dry run (${exec_mode} mode)"
        echo "--- ORCA input: ${inp_file} ---"
        cat "$inp_file"
        if [[ $exec_mode == slurm ]]; then
            local slurm_file="${workdir}/${tag}.slurm"
            write_slurm "$slurm_file" "$tag"
            echo ""
            echo "--- SLURM script: ${slurm_file} ---"
            cat "$slurm_file"
        else
            echo ""
            echo "--- Would run: ${orca_bin} ${inp_file} > ${log_file} ---"
        fi
        echo ""
        return
    fi

    # ---- live execution ------------------------------------------------
    if [[ $exec_mode == slurm ]]; then
        local slurm_file="${workdir}/${tag}.slurm"
        write_slurm "$slurm_file" "$tag"
        log "[${tag}] submitting to SLURM (partition=${partition})"
        (cd "$workdir" && sbatch "$(basename "$slurm_file")")
    else
        log "[${tag}] running ORCA locally (${cpus} cores)"
        "$orca_bin" "$inp_file" > "$log_file" 2>&1
        if grep -q "HURRAY\|THE OPTIMIZATION HAS CONVERGED" "$log_file" 2>/dev/null; then
            log "[${tag}] optimisation converged"
        else
            warn "[${tag}] convergence string not found — check ${log_file}"
        fi
    fi
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    local disp
    disp=$(disp_keyword)
    cat >&2 <<EOF
=============================================================
 Stage 1: Gas-Phase Geometry Optimisation
-------------------------------------------------------------
 Mode        : ${exec_mode}
 ORCA binary : ${orca_bin}
 Method      : ${method}${disp}
 Basis       : ${basis}
 DEFGRID     : ${grid}
 SCF MaxIter : ${max_iter}
 Cores       : ${cpus}
 MaxCore     : ${mem_mb} MB/core
 Dry run     : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    parse_cli "$@"

    # resolve execution mode and ORCA binary
    exec_mode=$(detect_mode "$force_local")
    orca_bin=$(resolve_orca_bin "$orca_bin_flag")

    # in non-dry-run mode, ORCA must be findable
    if ! $dry_run && [[ -z $orca_bin ]]; then
        die "ORCA binary not found. Set ORCA_BIN, use --orca-bin, or add orca to PATH."
    fi
    # even in dry-run, warn if missing
    if [[ -z $orca_bin ]]; then
        warn "ORCA binary not found — dry-run will proceed without it"
        orca_bin="orca"   # placeholder for display
    fi

    print_banner

    # process molecules from --list or single positional argument
    if [[ -n $list_file ]]; then
        while IFS= read -r entry || [[ -n $entry ]]; do
            [[ -z $entry || $entry == \#* || $entry == \;* ]] && continue
            local tag xyz_path
            if [[ $entry == *.xyz ]]; then
                xyz_path=$entry
                tag=$(basename "$entry" .xyz)
            else
                xyz_path="${XYZ_DIR}/${entry}.xyz"
                tag=$entry
            fi
            if [[ ! -f $xyz_path ]]; then
                warn "'${xyz_path}' not found — skipping"
                continue
            fi
            process_molecule "$xyz_path" "$tag"
        done < "$list_file"
    else
        local tag
        tag=$(basename "$xyz_arg" .xyz)
        process_molecule "$xyz_arg" "$tag"
    fi

    log "Stage 1 complete."
}

main "$@"
