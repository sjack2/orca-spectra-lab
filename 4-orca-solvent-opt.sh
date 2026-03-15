#!/usr/bin/env bash
# ============================================================================
# 4-orca-solvent-opt.sh -- Solvent-phase geometry optimization (Stage 4)
# ============================================================================
#
# OVERVIEW
#   For every conformer produced by Stage 3 (or 3b), this script:
#     1. Reads charge/multiplicity from the original XYZ in pre_xyz/,
#     2. Writes an ORCA input with implicit solvent (SMD or pure CPCM),
#     3. Either runs ORCA directly (--local) or submits a SLURM job array.
#
#   Each conformer is optimized independently. In local mode the
#   conformers are processed sequentially (this can take a while for
#   large ensembles -- use --dry-run first to check how many jobs there are).
#
# Usage:
#   4-orca-solvent-opt.sh TAG
#   4-orca-solvent-opt.sh --local --cpus 4 --solvent water ephedrine
#   4-orca-solvent-opt.sh --list mols.txt --dry-run
#
# Flags:
#   -m | --method NAME         DFT functional                  [B3LYP]
#   -b | --basis NAME          Basis set                       [def2-TZVP def2/J]
#        --disp KW             Dispersion: auto|none|D3BJ|D4   [auto]
#        --solvent NAME        Solvent name for implicit model   [water]
#        --solvent-model NAME  Implicit solvent model: smd|cpcm  [smd]
#        --max-iter N          SCF iteration limit              [150]
#   -c | --cpus N              CPU cores (%pal + SLURM)        [4]
#   -g | --grid N              DEFGRID level (1-3)             [3]
#        --mem-per-cpu MB      Memory per core in MB            [2048]
#        --orca-bin PATH       Path to ORCA binary              [auto]
#        --list FILE           File of molecule TAGs
#        --local               Run ORCA directly (no SLURM)
#        --dry-run             Write inputs but do not run
#   -h | --help                Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --partition NAME      SLURM partition                  [general]
#        --time HH:MM:SS       Wall-clock limit                 [06:00:00]
#        --max-running N       Max concurrent array tasks       [10]
#
# Cluster configuration (cluster.cfg):
#   Create a file named cluster.cfg in the same directory as this script
#   to set site-specific defaults without using flags every time.
#   See cluster.cfg.example for the full list of supported variables.
#   Variables set by cluster.cfg are overridden by explicit CLI flags.
#
# XYZ charge/multiplicity:
#   Parsed from pre_xyz/<TAG>.xyz line 2. Supports both formats:
#     charge=0 mult=1       (key=value)
#     0 1                   (bare integers)
#
# Directory layout:
#   <TAG>/
#   |-- 02_conf_search/split_xyz/     <- input conformers (from Stage 3/3b)
#   |-- 03_solvent_opt/
#   |   |-- <TAG>_conf_list.txt       conformer working-dir list (array input)
#   |   |-- <TAG>_array.slurm         SLURM array job script
#   |   |-- <TAG>_001/               one directory per conformer
#   |   |   |-- <TAG>_001.inp
#   |   |   |-- <TAG>_001.log
#   |   |   -- <TAG>_001.xyz        optimized geometry (ORCA output)
#   |   -- <TAG>_002/
#   |       -- ...
#   -- pre_xyz/<TAG>.xyz            must exist for charge/mult
#
# Examples:
#   # Local: optimise all conformers of ephedrine in water
#   4-orca-solvent-opt.sh --local --cpus 4 --solvent water ephedrine
#
#   # HPC: submit all conformers as a throttled job array (10 at a time)
#   4-orca-solvent-opt.sh --solvent acetonitrile --partition gpu aspirin
#
#   # HPC: raise concurrency limit for a small ensemble
#   4-orca-solvent-opt.sh --max-running 4 ephedrine
#
#   # Dry run: see how many conformers would be processed
#   4-orca-solvent-opt.sh --dry-run --local --list molecules.txt
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_METHOD="B3LYP"
DEFAULT_BASIS="def2-TZVP def2/J"
DEFAULT_DISP="auto"
DEFAULT_SOLVENT="water"
DEFAULT_SOLVENT_MODEL="smd"
DEFAULT_MAX_ITER=150
DEFAULT_CPUS=4
DEFAULT_GRID=3
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="general"
DEFAULT_WALL="06:00:00"
DEFAULT_MAX_RUNNING=10

XYZ_DIR="pre_xyz"
CONF_SUBDIR="02_conf_search/split_xyz"

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

# ============================================================================
# CLUSTER CONFIG
# ============================================================================
source_cluster_cfg() {
    local script_dir cfg
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cfg="${script_dir}/cluster.cfg"
    if [[ -f $cfg ]]; then
        log "Loaded cluster config: ${cfg}"
        # shellcheck source=/dev/null
        source "$cfg"
    fi
}

# ORCA binary resolution (same as Script 1)
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
    solvent=$DEFAULT_SOLVENT
    solvent_model=$DEFAULT_SOLVENT_MODEL
    max_iter=$DEFAULT_MAX_ITER
    cpus=$DEFAULT_CPUS
    grid=$DEFAULT_GRID
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=${CLUSTER_PARTITION:-$DEFAULT_PARTITION}
    wall=${CLUSTER_WALL:-$DEFAULT_WALL}
    max_running=${CLUSTER_MAX_RUNNING:-$DEFAULT_MAX_RUNNING}
    dry_run=false
    force_local=false
    list_file=""
    single_tag=""
    orca_bin_flag=""
    ompi_dir_flag=""

    local opts
    opts=$(getopt -o hb:m:c:g: \
        --long help,basis:,method:,disp:,solvent:,solvent-model:,max-iter:,cpus:,grid:,\
mem-per-cpu:,partition:,time:,max-running:,list:,orca-bin:,openmpi-dir:,local,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)      method=$2;         shift 2 ;;
            -b|--basis)       basis=$2;          shift 2 ;;
            --disp)           disp_mode=$2;      shift 2 ;;
            --solvent)        solvent=$2;         shift 2 ;;
            --solvent-model)  solvent_model=$2;  shift 2 ;;
            --max-iter)       max_iter=$2;       shift 2 ;;
            -c|--cpus)        cpus=$2;           shift 2 ;;
            -g|--grid)        grid=$2;           shift 2 ;;
            --mem-per-cpu)    mem_mb=$2;         shift 2 ;;
            --partition)      partition=$2;       shift 2 ;;
            --time)           wall=$2;           shift 2 ;;
            --max-running)    max_running=$2;    shift 2 ;;
            --list)           list_file=$2;      shift 2 ;;
            --orca-bin)       orca_bin_flag=$2;  shift 2 ;;
            --openmpi-dir)    ompi_dir_flag=$2;  shift 2 ;;
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
    [[ $solvent_model == smd || $solvent_model == cpcm ]] \
        || die "--solvent-model must be 'smd' or 'cpcm'"
}

# ============================================================================
# ORCA INPUT WRITER
# ============================================================================
write_orca_input() {
    local cid=$1 xyz_file=$2 inp_file=$3
    local disp cpcm_kw
    disp=$(disp_keyword)

    # SMD uses the bare CPCM keyword + %cpcm block; pure CPCM uses inline solvent name
    if [[ $solvent_model == smd ]]; then
        cpcm_kw="CPCM"
    else
        cpcm_kw="CPCM(${solvent})"
    fi

    cat >"$inp_file" <<EOF
! ${method} ${basis} Opt TightSCF RIJCOSX DEFGRID${grid}${disp} ${cpcm_kw}

%scf
    MaxIter ${max_iter}
end

%pal nprocs ${cpus} end
%maxcore ${mem_mb}

EOF

    if [[ $solvent_model == smd ]]; then
        cat >>"$inp_file" <<EOF
%cpcm
    SMD true
    SMDsolvent "${solvent}"
end

EOF
    fi

    cat >>"$inp_file" <<EOF
* xyz ${charge} ${mult}
$(tail -n +3 "$xyz_file")
*
EOF
}

# ============================================================================
# SLURM ARRAY JOB WRITER
# ============================================================================
write_array_slurm() {
    local tag=$1 conf_list=$2 slurm_file=$3 n=$4
    local abs_list out_dir orca_dir ompi_dir
    abs_list=$(cd "$(dirname "$conf_list")" && pwd)/$(basename "$conf_list")
    out_dir=$(dirname "$abs_list")
    orca_dir=$(dirname "$orca_bin")
    ompi_dir=${ompi_dir_flag:-${OMPI_DIR:-/shares/chem_hlw/orca/openmpi-4.1.6}}

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=solv_${tag}
#SBATCH --partition=${partition}
#SBATCH --array=1-${n}%${max_running}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${cpus}
#SBATCH --mem-per-cpu=${mem_mb}
#SBATCH --time=${wall}
#SBATCH --output=${out_dir}/slurm-%A_%a.out
#SBATCH --error=${out_dir}/slurm-%A_%a.err

# ---- ORCA / OpenMPI environment ----
export PATH="${ompi_dir}/bin:${orca_dir}:\$PATH"
export LD_LIBRARY_PATH="${ompi_dir}/lib:${orca_dir}:\$LD_LIBRARY_PATH"
export OPAL_PREFIX="${ompi_dir}"
export OMPI_MCA_btl="^openib"

WORKDIR=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${abs_list}")
CID=\$(basename "\$WORKDIR")
cd "\$WORKDIR"
"${orca_bin}" "\${CID}.inp" > "\${CID}.log"
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# PROCESS ONE CONFORMER  (local mode only)
# ============================================================================
process_conformer() {
    local tag=$1 cid=$2 xyz_file=$3

    local dir="${tag}/03_solvent_opt/${cid}"
    mkdir -p "$dir"

    local inp_file="${dir}/${cid}.inp"
    local log_file="${dir}/${cid}.log"

    write_orca_input "$cid" "$xyz_file" "$inp_file"

    if $dry_run; then
        log "  [${cid}] dry run -- input written to ${inp_file}"
        return
    fi

    log "  [${cid}] running ORCA (solvent=${solvent})"
    "$orca_bin" "$inp_file" > "$log_file" 2>&1
    if grep -q "HURRAY\|THE OPTIMIZATION HAS CONVERGED" "$log_file" 2>/dev/null; then
        log "  [${cid}] converged"
    else
        warn "  [${cid}] convergence not confirmed -- check ${log_file}"
    fi
}

# ============================================================================
# PROCESS ONE MOLECULE (all its conformers)
# ============================================================================
process_tag() {
    local tag=$1
    local xyz_dir="${tag}/${CONF_SUBDIR}"

    if [[ ! -d $xyz_dir ]]; then
        warn "[${tag}] conformer directory not found: ${xyz_dir}"
        return
    fi

    # read charge/multiplicity from the original pre_xyz file
    local pre_xyz="${XYZ_DIR}/${tag}.xyz"
    charge=0; mult=1
    if [[ -f $pre_xyz ]]; then
        read_xyz_header "$pre_xyz"
    else
        warn "[${tag}] ${pre_xyz} not found -- using charge=0 mult=1"
    fi

    # collect conformer XYZ files
    shopt -s nullglob
    local xyz_files=("${xyz_dir}/${tag}_"*.xyz)
    shopt -u nullglob

    local n=${#xyz_files[@]}
    if (( n == 0 )); then
        warn "[${tag}] no conformer XYZ files in ${xyz_dir}"
        return
    fi

    log "[${tag}] found ${n} conformers (charge=${charge}, mult=${mult})"

    if [[ $exec_mode == slurm ]]; then
        # Phase 1: write all ORCA inputs, collect absolute working dirs
        local out_dir="${tag}/03_solvent_opt"
        mkdir -p "$out_dir"
        local conf_dirs=()
        for xyz in "${xyz_files[@]}"; do
            local cid
            cid=$(basename "$xyz" .xyz)
            local dir="${tag}/03_solvent_opt/${cid}"
            mkdir -p "$dir"
            write_orca_input "$cid" "$xyz" "${dir}/${cid}.inp"
            if $dry_run; then
                log "  [${cid}] dry run -- input written to ${dir}/${cid}.inp"
            else
                conf_dirs+=("$(cd "$dir" && pwd)")
            fi
        done

        $dry_run && return

        # Phase 2: write conformer list, submit single throttled array job
        local conf_list="${out_dir}/${tag}_conf_list.txt"
        printf '%s\n' "${conf_dirs[@]}" > "$conf_list"

        local array_slurm="${out_dir}/${tag}_array.slurm"
        write_array_slurm "$tag" "$conf_list" "$array_slurm" "$n"
        local jid
        jid=$(sbatch --parsable "$array_slurm")
        log "[${tag}] submitted array job ${jid} (${n} conformers, max ${max_running} concurrent)"
    else
        for xyz in "${xyz_files[@]}"; do
            local cid
            cid=$(basename "$xyz" .xyz)
            process_conformer "$tag" "$cid" "$xyz"
        done
    fi
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    local disp smodel_label
    disp=$(disp_keyword)
    [[ $solvent_model == smd ]] && smodel_label="SMD via CPCM" || smodel_label="CPCM"
    cat >&2 <<EOF
=============================================================
 Stage 4: Solvent-Phase Geometry Optimization
-------------------------------------------------------------
 Mode        : ${exec_mode}
 ORCA binary : ${orca_bin}
 Method      : ${method}${disp}
 Basis       : ${basis}
 Solvent     : ${solvent} (${smodel_label})
 DEFGRID     : ${grid}
 SCF MaxIter : ${max_iter}
 Cores       : ${cpus}
 MaxCore     : ${mem_mb} MB/core
 Max running : ${max_running} (SLURM array throttle)
 Dry run     : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"

    exec_mode=$(detect_mode "$force_local")
    orca_bin=$(resolve_orca_bin "$orca_bin_flag")

    if ! $dry_run && [[ -z $orca_bin ]]; then
        die "ORCA binary not found. Set ORCA_BIN in cluster.cfg, use --orca-bin, or add orca to PATH."
    fi
    if ! $dry_run && [[ -n $orca_bin && ! -x $orca_bin ]]; then
        die "ORCA binary '${orca_bin}' does not exist or is not executable. Check ORCA_BIN in cluster.cfg or --orca-bin."
    fi
    if [[ -z $orca_bin ]]; then
        warn "ORCA binary not found -- dry-run will proceed without it"
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

    log "Stage 4 complete."
}

main "$@"
