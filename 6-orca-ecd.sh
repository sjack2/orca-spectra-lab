#!/usr/bin/env bash
# ============================================================================
# 6-orca-ecd.sh — TD-DFT excited-state calculations (Stage 6)
# ============================================================================
#
# OVERVIEW
#   For every conformer that passed the Boltzmann filter (Stage 5), this
#   script builds an ORCA TD-DFT input with implicit solvent and either
#   runs it locally or submits a SLURM job.
#
#   The resulting outputs contain UV-Vis oscillator strengths and ECD
#   rotatory strengths for all requested excited states.
#
# Usage:
#   6-orca-ecd.sh TAG
#   6-orca-ecd.sh --local --cpus 4 --method wB97X-D3 ephedrine
#   6-orca-ecd.sh --list mols.txt --dry-run
#
# Flags:
#   -m | --method NAME         DFT functional                  [wB97X-D3]
#   -b | --basis NAME          Basis set                       [def2-TZVP def2/J]
#        --disp KW             Dispersion: auto|none|D3BJ|D4   [auto]
#        --roots N             Number of excited states         [30]
#        --solvent NAME        SMD solvent keyword              [water]
#        --max-iter N          SCF iteration limit              [150]
#   -c | --cpus N              CPU cores (%pal + SLURM)        [4]
#   -g | --grid N              DEFGRID level (1–3)             [3]
#        --mem-per-cpu MB      Memory per core in MB            [2048]
#        --orca-bin PATH       Path to ORCA binary              [auto]
#        --list FILE           File of molecule TAGs
#        --local               Run ORCA directly (no SLURM)
#        --dry-run             Write inputs but do not run
#   -h | --help                Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --partition NAME      SLURM partition                  [circe]
#        --time HH:MM:SS       Wall-clock limit                 [06:00:00]
#
# Directory layout:
#   <TAG>/
#   ├── bw_results/<TAG>_bw_labels.dat     ← conformer list (from Stage 5)
#   ├── solvent_opt/<CID>/<CID>.xyz        ← geometries (from Stage 4)
#   └── ecd/
#       ├── <CID>/
#       │   ├── <CID>.inp                  ORCA TD-DFT input
#       │   ├── <CID>.log                  ORCA output
#       │   └── <CID>.slurm               (HPC mode only)
#       └── ...
#
# Examples:
#   # Local: compute ECD for all populated conformers of ephedrine
#   6-orca-ecd.sh --local --cpus 4 --roots 30 ephedrine
#
#   # HPC: submit with a range-separated hybrid
#   6-orca-ecd.sh --method CAM-B3LYP --partition gpu aspirin
#
#   # Dry run: inspect inputs
#   6-orca-ecd.sh --dry-run --local --list molecules.txt
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_METHOD="wB97X-D3"
DEFAULT_BASIS="def2-TZVP def2/J"
DEFAULT_DISP="auto"
DEFAULT_ROOTS=30
DEFAULT_SOLVENT="water"
DEFAULT_MAX_ITER=150
DEFAULT_CPUS=4
DEFAULT_GRID=3
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="circe"
DEFAULT_WALL="06:00:00"

XYZ_DIR="pre_xyz"

# ============================================================================
# HELPERS
# ============================================================================
die()  { printf 'Error: %s\n' "$*" >&2; exit 1; }
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
warn() { printf '[%s] Warning: %s\n' "$(date '+%F %T')" "$*" >&2; }

require_file() { [[ -f $1 ]] || die "File '$1' not found"; }

show_help() {
    sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# \{0,1\}//'
    exit 0
}

resolve_orca_bin() {
    local flag_val=$1
    if [[ -n $flag_val ]]; then printf '%s' "$flag_val"; return; fi
    if [[ -n ${ORCA_BIN:-} ]]; then printf '%s' "$ORCA_BIN"; return; fi
    local found; found=$(command -v orca 2>/dev/null || true)
    if [[ -n $found ]]; then printf '%s' "$found"; return; fi
    printf ''
}

detect_mode() {
    local force_local=$1
    if $force_local; then printf 'local'
    elif command -v sbatch >/dev/null 2>&1; then printf 'slurm'
    else printf 'local'; fi
}

# ============================================================================
# CHARGE / MULTIPLICITY PARSER
# ============================================================================
read_xyz_header() {
    local xyz_file=$1
    local header_line
    header_line=$(sed -n '2p' "$xyz_file")

    if [[ $header_line =~ charge[[:space:]]*=[[:space:]]*([+-]?[0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
    else
        charge=0
    fi

    if [[ $header_line =~ mult[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
        mult=${BASH_REMATCH[1]}
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
disp_keyword() {
    local method_upper=${method^^}
    case $disp_mode in
        none|NONE) printf '' ;;
        D3BJ|d3bj) printf ' D3BJ' ;;
        D4|d4) printf ' D4' ;;
        auto|AUTO)
            if [[ $method =~ (-D[0-9]?|-D3BJ|-D3ZERO|-D4)($|[[:space:]]) ]]; then
                printf ''; return; fi
            case $method_upper in
                WB97X-D|WB97X-D3|WB97X-D4|WB97X-V|WB97XD|WB97M-V|B97-D|B97-D3)
                    printf ''; return ;; esac
            printf ' D3BJ' ;;
        *) die "--disp must be auto, none, D3BJ, or D4" ;;
    esac
}

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    method=$DEFAULT_METHOD
    basis=$DEFAULT_BASIS
    disp_mode=$DEFAULT_DISP
    roots=$DEFAULT_ROOTS
    solvent=$DEFAULT_SOLVENT
    max_iter=$DEFAULT_MAX_ITER
    cpus=$DEFAULT_CPUS
    grid=$DEFAULT_GRID
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=$DEFAULT_PARTITION
    wall=$DEFAULT_WALL
    dry_run=false
    force_local=false
    list_file=""
    single_tag=""
    orca_bin_flag=""

    local opts
    opts=$(getopt -o hb:m:c:g: \
        --long help,basis:,method:,disp:,roots:,solvent:,max-iter:,cpus:,grid:,\
mem-per-cpu:,partition:,time:,list:,orca-bin:,local,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)      method=$2;         shift 2 ;;
            -b|--basis)       basis=$2;          shift 2 ;;
            --disp)           disp_mode=$2;      shift 2 ;;
            --roots)          roots=$2;          shift 2 ;;
            --solvent)        solvent=$2;         shift 2 ;;
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

    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional args not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide exactly one molecule TAG"
        single_tag=$1
    fi

    (( grid < 1 )) && grid=1 || true
    (( grid > 3 )) && grid=3 || true
}

# ============================================================================
# ORCA INPUT WRITER
# ============================================================================
write_orca_input() {
    local cid=$1 xyz_file=$2 inp_file=$3
    local disp
    disp=$(disp_keyword)

    cat >"$inp_file" <<EOF
! ${method} ${basis} TightSCF RIJCOSX DEFGRID${grid}${disp} CPCM

%scf
    MaxIter ${max_iter}
end

%pal nprocs ${cpus} end
%maxcore ${mem_mb}

%tddft
    nroots ${roots}
    TDA false
end

%cpcm
    SMD true
    SMDsolvent "${solvent}"
end

* xyz ${charge} ${mult}
$(tail -n +3 "$xyz_file")
*
EOF
}

# ============================================================================
# SLURM WRITER
# ============================================================================
write_slurm() {
    local cid=$1 inp_file=$2 slurm_file=$3
    local orca_dir
    orca_dir=$(dirname "$orca_bin")

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=ecd_${cid}
#SBATCH --partition=${partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${cpus}
#SBATCH --mem-per-cpu=${mem_mb}
#SBATCH --time=${wall}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

export PATH="${orca_dir}:\$PATH"
export LD_LIBRARY_PATH="${orca_dir}:\$LD_LIBRARY_PATH"

"${orca_bin}" "${cid}.inp" > "${cid}.log"
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# PROCESS ONE CONFORMER
# ============================================================================
process_conformer() {
    local tag=$1 cid=$2 xyz_file=$3

    local dir="${tag}/ecd/${cid}"
    mkdir -p "$dir"

    local inp_file="${dir}/${cid}.inp"
    local log_file="${dir}/${cid}.log"

    write_orca_input "$cid" "$xyz_file" "$inp_file"

    if $dry_run; then
        log "  [${cid}] dry run — input written to ${inp_file}"
        return
    fi

    if [[ $exec_mode == slurm ]]; then
        local slurm_file="${dir}/${cid}.slurm"
        write_slurm "$cid" "$inp_file" "$slurm_file"
        (cd "$dir" && sbatch "$(basename "$slurm_file")")
    else
        log "  [${cid}] running TD-DFT (${roots} roots, solvent=${solvent})"
        "$orca_bin" "$inp_file" > "$log_file" 2>&1
        if grep -q "ABSORPTION SPECTRUM" "$log_file" 2>/dev/null; then
            log "  [${cid}] TD-DFT completed successfully"
        else
            warn "  [${cid}] spectrum block not found — check ${log_file}"
        fi
    fi
}

# ============================================================================
# PROCESS ONE MOLECULE (all its Boltzmann-filtered conformers)
# ============================================================================
process_tag() {
    local tag=$1

    # find the label file from Stage 5
    local lab_file="${tag}/bw_results/${tag}_bw_labels.dat"
    if [[ ! -f $lab_file ]]; then
        warn "[${tag}] label file not found: ${lab_file}"
        return
    fi

    # read conformer IDs
    local cids=()
    mapfile -t cids < <(grep -v '^[[:space:]]*$' "$lab_file")
    if (( ${#cids[@]} == 0 )); then
        warn "[${tag}] no conformers in label file"
        return
    fi

    # read charge/multiplicity from the original XYZ
    charge=0; mult=1
    local pre_xyz="${XYZ_DIR}/${tag}.xyz"
    if [[ -f $pre_xyz ]]; then
        read_xyz_header "$pre_xyz"
    else
        warn "[${tag}] ${pre_xyz} not found — using charge=0 mult=1"
    fi

    log "[${tag}] ${#cids[@]} conformers to process (charge=${charge}, mult=${mult})"

    for cid in "${cids[@]}"; do
        cid=$(echo "$cid" | xargs)  # trim whitespace
        [[ -z $cid ]] && continue

        # look for the solvent-optimised geometry
        local xyz="${tag}/solvent_opt/${cid}/${cid}.xyz"
        if [[ ! -f $xyz ]]; then
            warn "  [${cid}] solvent-optimised XYZ not found — skipping"
            continue
        fi
        process_conformer "$tag" "$cid" "$xyz"
    done
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    local disp
    disp=$(disp_keyword)
    cat >&2 <<EOF
=============================================================
 Stage 6: TD-DFT Excited-State Calculation (UV-Vis + ECD)
-------------------------------------------------------------
 Mode        : ${exec_mode}
 ORCA binary : ${orca_bin}
 Method      : ${method}${disp}
 Basis       : ${basis}
 Roots       : ${roots}
 TDA         : false (full TD-DFT)
 Solvent     : ${solvent} (SMD via CPCM)
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

    exec_mode=$(detect_mode "$force_local")
    orca_bin=$(resolve_orca_bin "$orca_bin_flag")

    if ! $dry_run && [[ -z $orca_bin ]]; then
        die "ORCA binary not found. Set ORCA_BIN, use --orca-bin, or add orca to PATH."
    fi
    if [[ -z $orca_bin ]]; then
        warn "ORCA binary not found — dry-run will proceed without it"
        orca_bin="orca"
    fi

    print_banner

    if [[ -n $list_file ]]; then
        while IFS= read -r tag || [[ -n $tag ]]; do
            [[ -z $tag || $tag == \#* || $tag == \;* ]] && continue
            process_tag "$tag"
        done < "$list_file"
    else
        process_tag "$single_tag"
    fi

    log "Stage 6 complete."
}

main "$@"
