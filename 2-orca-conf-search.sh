#!/usr/bin/env bash
# ============================================================================
# 2-orca-conf-search.sh — Conformer enumeration with Open Babel Confab
# ============================================================================
#
# OVERVIEW
#   Stage 2 of the ORCA workflow.  For each molecule it:
#     1. Retrieves the optimised geometry from Stage 1 (or a user-supplied XYZ),
#     2. Converts that geometry to SDF format,
#     3. Runs Open Babel Confab to generate conformers within an energy window.
#
#   Confab operates by systematically rotating torsion angles and evaluating
#   each resulting geometry with a force field.  Structures above the energy
#   cutoff are discarded.  For rigid molecules Confab may return only a single
#   conformer; for flexible molecules it can return hundreds.
#
# Usage:
#   2-orca-conf-search.sh TAG
#   2-orca-conf-search.sh path/to/foo.xyz --conf 500
#   2-orca-conf-search.sh --list molecules.txt --ecut 4 --dry-run
#
# Flags:
#   --conf N       Maximum number of conformers              [1000]
#   --ecut E       Energy cutoff in kcal/mol                 [6]
#   --list FILE    Text file of TAGs or XYZ paths
#   --dry-run      Echo commands without executing them
#   -h | --help    Show this help and exit
#
# Directory layout (reads from Stage 1, writes to):
#   <TAG>/
#   ├── <TAG>_orca_opt/
#   │   └── <TAG>.xyz              ← input (from Stage 1)
#   └── orca_opt_conf/
#       ├── <TAG>.xyz              copy of optimised geometry
#       ├── <TAG>.sdf              SDF conversion
#       ├── <TAG>_combined.sdf     all conformers in one file
#       ├── split_sdf/             individual SDF files
#       └── split_xyz/             individual XYZ files  → Stage 3 input
#
# Examples:
#   2-orca-conf-search.sh --ecut 10 --conf 500 ephedrine
#   2-orca-conf-search.sh --list molecules.txt --dry-run
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_CONF_COUNT=1000
DEFAULT_ECUT_KCAL=6
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

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    conf_count=$DEFAULT_CONF_COUNT
    ecut_kcal=$DEFAULT_ECUT_KCAL
    dry_run=false
    list_file=""
    positional=()

    local opts
    opts=$(getopt -o h --long help,list:,conf:,ecut:,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            --list)    list_file=$2; shift 2 ;;
            --conf)    conf_count=$2; shift 2 ;;
            --ecut)    ecut_kcal=$2; shift 2 ;;
            --dry-run) dry_run=true; shift ;;
            -h|--help) show_help ;;
            --)        shift; break ;;
            *)         die "Unknown option '$1'" ;;
        esac
    done
    positional+=("$@")

    if [[ -n $list_file ]]; then
        [[ ${#positional[@]} -eq 0 ]] || die "Positional args invalid with --list"
        require_file "$list_file"
    else
        [[ ${#positional[@]} -gt 0 ]] || die "Provide at least one TAG or XYZ path"
    fi

    command -v obabel >/dev/null 2>&1 || die "obabel not found in PATH"
}

# ============================================================================
# PATH RESOLUTION
# ============================================================================
# Look for the optimised XYZ from Stage 1 in the expected location.
find_optimised_xyz() {
    local tag=$1
    local probe="${tag}/${tag}_orca_opt/${tag}.xyz"
    if [[ -f $probe ]]; then
        printf '%s' "$probe"
    else
        printf ''
    fi
}

# ============================================================================
# CORE: CONFAB CONFORMER GENERATION
# ============================================================================
process_tag() {
    local tag=$1 xyz_path=$2

    local out_dir="${tag}/orca_opt_conf"
    mkdir -p "$out_dir"

    # copy optimised geometry into the conformer directory
    local xyz_copy="${out_dir}/${tag}.xyz"
    cp -f "$xyz_path" "$xyz_copy"

    local sdf="${out_dir}/${tag}.sdf"
    local combined="${out_dir}/${tag}_combined.sdf"

    if $dry_run; then
        log "[${tag}] (dry run) obabel -ixyz ${xyz_copy} -osdf -O ${sdf}"
        log "[${tag}] (dry run) obabel ${sdf} -O ${combined} --confab --original --conf ${conf_count} --ecutoff ${ecut_kcal}"
        return
    fi

    # convert XYZ → SDF (Confab requires SDF input)
    log "[${tag}] converting to SDF"
    obabel -ixyz "$xyz_copy" -osdf -O "$sdf" 2>/dev/null
    if [[ ! -s $sdf ]]; then
        warn "[${tag}] SDF conversion produced empty file — skipping"
        return
    fi

    # run Confab
    log "[${tag}] running Confab (ecutoff=${ecut_kcal} kcal/mol, max=${conf_count})"
    obabel "$sdf" -O "$combined" \
        --confab --original --conf "$conf_count" --ecutoff "$ecut_kcal" 2>/dev/null

    if [[ ! -s $combined ]]; then
        warn "[${tag}] Confab produced no conformers"
        return
    fi

    # split combined SDF into individual files and convert to XYZ
    local split_sdf="${out_dir}/split_sdf"
    local split_xyz="${out_dir}/split_xyz"
    mkdir -p "$split_sdf" "$split_xyz"

    log "[${tag}] splitting conformers"
    obabel "$combined" -O "${split_sdf}/${tag}_".sdf -m 2>/dev/null

    shopt -s nullglob
    local sdf_files=("${split_sdf}/${tag}_"*.sdf)
    shopt -u nullglob

    local count=${#sdf_files[@]}
    if (( count == 0 )); then
        warn "[${tag}] no split SDF files generated"
        return
    fi

    for s in "${sdf_files[@]}"; do
        local base
        base=$(basename "${s%.sdf}")
        obabel "$s" -O "${split_xyz}/${base}.xyz" 2>/dev/null
    done

    log "[${tag}] generated ${count} conformers"
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 2: Conformer Generation (Confab)
-------------------------------------------------------------
 Energy cutoff : ${ecut_kcal} kcal/mol
 Max conformers: ${conf_count}
 Dry run       : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    parse_cli "$@"
    print_banner

    local entries=()
    if [[ -n $list_file ]]; then
        mapfile -t entries < <(grep -vE '^[[:space:]]*(#|;|$)' "$list_file")
    else
        entries=("${positional[@]}")
    fi

    for entry in "${entries[@]}"; do
        local tag xyz_path
        if [[ $entry == *.xyz ]]; then
            tag=$(basename "$entry" .xyz)
            xyz_path=$entry
            require_file "$xyz_path"
        else
            tag=$entry
            xyz_path=$(find_optimised_xyz "$tag")
            if [[ -z $xyz_path ]]; then
                warn "[${tag}] optimised XYZ not found in ${tag}/${tag}_orca_opt/ — skipping"
                continue
            fi
        fi
        process_tag "$tag" "$xyz_path"
    done

    log "Stage 2 complete."
}

main "$@"
