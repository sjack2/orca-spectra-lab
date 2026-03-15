# ORCA Computational Spectra Lab

A modular Bash/Python workflow for computing UV-Vis absorption, electronic circular dichroism (ECD), and vibrational circular dichroism (VCD) spectra using ORCA density functional theory. Designed for teaching and research, and runs identically on a local workstation or a SLURM-managed HPC cluster.

> **Accompanying article:**  _A Computational Spectroscopy Module for Teaching DFT-Based Spectral Workflows in Graduate Chemistry_ (in preparation for the _Journal of Chemical Education_).

---

## Features

- **Six-stage pipeline** from raw XYZ geometry to publication-quality broadened spectra, with parallel Stage 6 branches for electronic (UV-Vis/ECD) and vibrational (IR/VCD) spectroscopy.
- **Dual execution mode:** every ORCA-calling script auto-detects whether SLURM is available. Pass `--local` to force direct execution; omit it on a cluster and jobs are submitted via `sbatch`.
- **Two conformer search paths:** Open Babel Confab (fast, force-field-based) or CREST (GFN2-xTB metadynamics, more thorough).
- **Boltzmann-weighted spectral averaging** with configurable temperature and population threshold.
- **Physically correct broadening:** Gaussian convolution in energy space (eV) for electronic spectra and in wavenumber space (cm-1) for vibrational spectra, with Jacobian correction for ECD. Produces both PNG/PDF plots and CSV data tables.
- **`--dry-run` on every script** to inspect generated ORCA inputs without running any calculations.

---

## Pipeline Overview

```
Stage 1                    Gas-phase geometry optimization
  |                        1-orca-init-opt.sh
  v
Stage 2                    Conformer enumeration
  |                        2-orca-conf-search.sh  (Confab)
  |                   or   2b-crest-conf-search.sh (CREST)
  v
Stage 3                    Split conformers into individual XYZ files
  |                        3-orca-conf-split.sh   (Confab output)
  |                   or   3b-crest-conf-split.sh  (CREST output)
  v
Stage 4                    Solvent-phase re-optimization (SMD/CPCM)
  |                        4-orca-solvent-opt.sh
  v
Stage 5                    Boltzmann weighting & filtering
  |                        5-orca-boltzmann-weight.sh
  |
  |------------------------------------------+
  v                                          v
Stage 6-ecd                              Stage 6-vcd
  TD-DFT excited-state calculations        Analytic frequency calculations
  6-orca-ecd.sh                            6-orca-vcd.sh
  |                                          |
  v                                          v
Plot                                     Plot
  Broadened UV-Vis & ECD spectra           Broadened IR & VCD spectra
  or_ecd_uvvis_tools.py                    or_vcd_ir_tools.py
```

Stages 1-5 are shared. After Boltzmann filtering, the pipeline branches: run `6-orca-ecd.sh` for electronic spectra, `6-orca-vcd.sh` for vibrational spectra, or both.

---

## Quick Start

### 1. Install prerequisites

See [INSTALL.md](INSTALL.md) for detailed instructions. In brief, you need:

| Software | Version | Purpose |
|----------|---------|---------|
| ORCA | 6.0+ | Quantum chemistry engine |
| OpenMPI | >= 4.0 | Parallel execution for ORCA |
| Open Babel | >= 3.0 | File conversion & Confab conformer search |
| Python 3 | >= 3.8 | Plotting tools |
| CREST | >= 2.12 | _(optional)_ GFN2-xTB conformer search |

### 2. Clone the repository

```bash
git clone https://github.com/sjack2/orca-spectra-lab.git
cd orca-spectra-lab
```

### 3. Check Python and install dependencies

First verify that Python 3.8 or later is available:

```bash
python3 --version
```

If the command is not found or the version is below 3.8, install Python before continuing. On an HPC cluster, the easiest approach is [Miniconda](https://docs.conda.io/en/latest/miniconda.html):

```bash
# Download and install Miniconda (Linux, x86-64)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow the prompts, then open a new shell or run: source ~/.bashrc
```

Once Python is available, install the plotting dependencies:

```bash
pip install -r requirements.txt
```

### 4. Make scripts executable

```bash
chmod +x *.sh
```

### 5. Configure for your cluster (HPC users)

Copy the provided template and fill in the paths for your site. This file is read automatically by every ORCA-calling script, so you will not need to pass `--orca-bin` or `--partition` on every invocation.

```bash
cp cluster.cfg.example cluster.cfg
```

Then edit `cluster.cfg` with your site's values, for example:

```bash
ORCA_BIN=/shares/chem_hlw/orca/orca_6_0_1_linux_x86-64_shared_openmpi416_avx2/orca
OMPI_DIR=/shares/chem_hlw/orca/openmpi-4.1.6
CLUSTER_PARTITION=general
CLUSTER_WALL=06:00:00
```

`cluster.cfg` is gitignored, so site-specific paths are never committed to the repository. The template `cluster.cfg.example` documents all supported variables.

### 6. Set up your molecule

Place a starting XYZ file in `pre_xyz/`. Line 2 must encode the charge and multiplicity:

```
14
charge=0 mult=1
C    0.000000    0.000000    0.000000
...
```

Both `charge=0 mult=1` (key=value) and `0 1` (bare integers) are accepted. If neither is found, defaults are charge=0, mult=1.

### 7. Run the pipeline (local workstation example)

```bash
# Stage 1: Optimize geometry in vacuum
./1-orca-init-opt.sh --local --cpus 4 pna

# Stage 2: Generate conformers (Confab path)
./2-orca-conf-search.sh --ecut 5 --conf 500 pna
# Stage 3: Split conformers into individual XYZ files
./3-orca-conf-split.sh pna

# Stage 4: Re-optimize each conformer in solvent
./4-orca-solvent-opt.sh --local --cpus 4 --solvent water pna

# Stage 5: Boltzmann filter
./5-orca-boltzmann-weight.sh pna

# --- Electronic spectra branch ---
# Stage 6-ecd: TD-DFT on populated conformers
./6-orca-ecd.sh --local --cpus 4 --method wB97X-D3 --roots 30 pna

# Plot UV-Vis & ECD spectra
python3 or_ecd_uvvis_tools.py pna/05_ecd \
    --bw pna/04_boltzmann/pna_energies.dat \
    --outdir pna --xlim 200 500

# --- Vibrational spectra branch ---
# Stage 6-vcd: Analytic frequencies on populated conformers
./6-orca-vcd.sh --local --cpus 4 --method B3LYP --solvent water pna

# Plot IR & VCD spectra
python3 or_vcd_ir_tools.py pna/05_vcd \
    --bw pna/04_boltzmann/pna_energies.dat \
    --outdir pna --xlim 800 3500
```

### 8. Run the pipeline (HPC/SLURM example)

On a cluster with SLURM, simply omit `--local`. The scripts detect `sbatch` automatically and submit jobs. With `cluster.cfg` configured, no additional flags are needed:

```bash
./1-orca-init-opt.sh --cpus 8 pna
# submits SLURM job to the partition set in cluster.cfg; monitor with squeue
```

---

## Directory Structure

Every script follows the same directory convention. For a molecule tagged `aspirin`:

```
aspirin/
|-- 01_gas_opt/                Stage 1 - gas-phase optimization
|   |-- aspirin.inp
|   |-- aspirin.log
|   `-- aspirin.xyz            optimized geometry
|-- 02_conf_search/            Stage 2 - conformer enumeration (Confab or CREST)
|   |-- aspirin.xyz            copy of Stage 1 geometry
|   |-- aspirin_combined.sdf   all conformers (Confab output)
|   |-- split_sdf/             individual SDF files (Stage 3 writes these)
|   `-- split_xyz/             individual conformer XYZ files (Stage 3 writes these)
|       |-- aspirin_1.xyz
|       |-- aspirin_2.xyz
|       `-- ...
|-- 03_solvent_opt/            Stage 4 - solvent-phase re-optimization
|   |-- aspirin_1/
|   |   |-- aspirin_1.inp
|   |   |-- aspirin_1.log
|   |   `-- aspirin_1.xyz      optimized geometry per conformer
|   `-- aspirin_2/
|       `-- ...
|-- 04_boltzmann/              Stage 5 - Boltzmann weighting
|   |-- aspirin_energies.dat   full table (CID, E, dE, p)
|   `-- aspirin_bw_labels.dat  conformer IDs above threshold
|-- 05_ecd/                    Stage 6-ecd - TD-DFT excited states
|   |-- aspirin_1/
|   |   |-- aspirin_1.inp
|   |   `-- aspirin_1.log
|   `-- aspirin_2/
|       `-- ...
|-- 05_vcd/                    Stage 6-vcd - analytic frequencies
|   |-- aspirin_1/
|   |   |-- aspirin_1.inp
|   |   `-- aspirin_1.log
|   `-- aspirin_2/
|       `-- ...
`-- 06_spectra/                Plotting outputs
    |-- aspirin_uvvis.png
    |-- aspirin_ecd.png
    |-- aspirin_ir.png
    |-- aspirin_vcd.png
    `-- *.csv / *.pdf
```

**Note:** Stage 2 only generates the combined SDF/XYZ ensemble. Stage 3 (or 3b) is always required to split that output into the per-conformer XYZ files consumed by Stage 4.

Starting geometries live in `pre_xyz/`:

```
pre_xyz/
|-- aspirin.xyz
|-- pna.xyz
|-- ephedrine.xyz
`-- methyloxirane.xyz
```

---

## CLI Reference

All scripts accept `--help` for full usage. Flags shown with `[default]`.

### 1-orca-init-opt.sh -- Gas-Phase Geometry Optimization

```
1-orca-init-opt.sh [OPTIONS] TAG [TAG ...]
1-orca-init-opt.sh [OPTIONS] --list FILE
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | DFT functional | `B3LYP` |
| `--basis NAME` | `-b` | Basis set(s) | `def2-SVP def2/J` |
| `--disp KW` | | Dispersion: `auto`, `none`, `D3BJ`, `D4` | `auto` |
| `--max-iter N` | | SCF iteration limit | `300` |
| `--cpus N` | `-c` | CPU cores | `4` |
| `--grid N` | `-g` | DEFGRID level (1-3) | `3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--orca-bin PATH` | | Path to ORCA binary | _auto-detected_ |
| `--list FILE` | | Text file of TAGs (one per line) | |
| `--local` | | Run ORCA directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Output:** `<TAG>/01_gas_opt/<TAG>.xyz` (optimized geometry)

### 2-orca-conf-search.sh -- Confab Conformer Search

```
2-orca-conf-search.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--conf N` | Maximum conformers | `1000` |
| `--ecut E` | Energy cutoff (kcal/mol) | `6` |
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo commands without running | |

**Output:** `<TAG>/02_conf_search/<TAG>_combined.sdf`

### 2b-crest-conf-search.sh -- CREST Conformer Search

```
2b-crest-conf-search.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--cpus N` | `-c` | CPU threads for CREST | `4` |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run CREST directly (no SLURM) | |
| `--dry-run` | | Echo actions without running | |
| `--pre-xyz` | | Also search `pre_xyz/` for input | |
| `--mem MB` | | SLURM memory _(HPC only)_ | `4096` |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall time _(HPC only)_ | `06:00:00` |

**Output:** `<TAG>/02_conf_search/crest_conformers.xyz`

### 3-orca-conf-split.sh -- Split Confab Output

```
3-orca-conf-split.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo actions without running | |

**Input:** `<TAG>/02_conf_search/<TAG>_combined.sdf`
**Output:** `<TAG>/02_conf_search/split_xyz/<TAG>_N.xyz`

### 3b-crest-conf-split.sh -- Split CREST Output

```
3b-crest-conf-split.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo actions without running | |

**Input:** `<TAG>/02_conf_search/crest_conformers.xyz`
**Output:** `<TAG>/02_conf_search/split_xyz/<TAG>_NNN.xyz` (zero-padded)

### 4-orca-solvent-opt.sh -- Solvent-Phase Re-optimization

```
4-orca-solvent-opt.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | DFT functional | `B3LYP` |
| `--basis NAME` | `-b` | Basis set(s) | `def2-TZVP def2/J` |
| `--disp KW` | | Dispersion: `auto`, `none`, `D3BJ`, `D4` | `auto` |
| `--solvent NAME` | | SMD solvent keyword | `water` |
| `--max-iter N` | | SCF iteration limit | `150` |
| `--cpus N` | `-c` | CPU cores | `4` |
| `--grid N` | `-g` | DEFGRID level (1-3) | `3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--orca-bin PATH` | | Path to ORCA binary | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run ORCA directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Input:** `<TAG>/02_conf_search/split_xyz/` (from Stage 3/3b)
**Output:** `<TAG>/03_solvent_opt/<CID>/<CID>.xyz` (optimized geometry per conformer)

### 5-orca-boltzmann-weight.sh -- Boltzmann Weighting & Filtering

```
5-orca-boltzmann-weight.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--temp K` | Temperature in Kelvin | `298.15` |
| `--p-cut VAL` | Probability cutoff (e.g., 0.01 = 1%) | `0.01` |
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Show what would be computed | |

**Input:** `<TAG>/03_solvent_opt/<CID>/<CID>.log` (ORCA outputs from Stage 4)
**Output:**
- `<TAG>/04_boltzmann/<TAG>_energies.dat` - full table: conformer ID, energy (Hartree), relative energy (kcal/mol), Boltzmann probability
- `<TAG>/04_boltzmann/<TAG>_bw_labels.dat` - conformer IDs above the probability cutoff

### 6-orca-ecd.sh -- TD-DFT Excited-State Calculations

```
6-orca-ecd.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | DFT functional | `wB97X-D3` |
| `--basis NAME` | `-b` | Basis set(s) | `def2-TZVP def2/J` |
| `--disp KW` | | Dispersion: `auto`, `none`, `D3BJ`, `D4` | `auto` |
| `--roots N` | | Number of excited states | `30` |
| `--solvent NAME` | | SMD solvent keyword | `water` |
| `--max-iter N` | | SCF iteration limit | `150` |
| `--cpus N` | `-c` | CPU cores | `4` |
| `--grid N` | `-g` | DEFGRID level (1-3) | `3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--orca-bin PATH` | | Path to ORCA binary | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run ORCA directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Input:** `<TAG>/04_boltzmann/<TAG>_bw_labels.dat` (from Stage 5) + `<TAG>/03_solvent_opt/<CID>/<CID>.xyz` (from Stage 4)
**Output:** `<TAG>/05_ecd/<CID>/<CID>.log` (ORCA TD-DFT output with absorption and CD spectrum blocks)

### 6-orca-vcd.sh -- Analytic Frequency / IR + VCD Calculations

```
6-orca-vcd.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | DFT functional | `B3LYP` |
| `--basis NAME` | `-b` | Basis set(s) | `def2-TZVP def2/J` |
| `--disp KW` | | Dispersion: `auto`, `none`, `D3BJ`, `D4` | `auto` |
| `--solvent NAME` | | SMD solvent keyword | `water` |
| `--max-iter N` | | SCF iteration limit | `150` |
| `--cpus N` | `-c` | CPU cores | `4` |
| `--grid N` | `-g` | DEFGRID level (1-3) | `3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--orca-bin PATH` | | Path to ORCA binary | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run ORCA directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Input:** `<TAG>/04_boltzmann/<TAG>_bw_labels.dat` (from Stage 5) + `<TAG>/03_solvent_opt/<CID>/<CID>.xyz` (from Stage 4)
**Output:** `<TAG>/05_vcd/<CID>/<CID>.log` (ORCA AnFreq output with IR and VCD spectrum blocks)

### or_ecd_uvvis_tools.py -- Electronic Spectral Broadening & Plotting

```
or_ecd_uvvis_tools.py [OPTIONS] LOGS [LOGS ...]
```

`LOGS` can be individual `.log` files, glob patterns, or directories (searched recursively).

| Flag | Description | Default |
|------|-------------|---------|
| `--bw PATH` | Boltzmann weight file (`_energies.dat` from Stage 5) | _equal weights_ |
| `--outdir TAG` | Molecule directory; outputs go to `<TAG>/06_spectra/<TAG>_*` | |
| `--prefix STR` | Explicit output prefix (overrides `--outdir`) | `spectra` |
| `--title STR` | Plot title (overrides auto-generated title) | _derived from prefix_ |
| `--uv_fwhm EV` | Gaussian FWHM for UV-Vis (eV) | `0.35` |
| `--ecd_fwhm EV` | Gaussian FWHM for ECD (eV) | `0.25` |
| `--xlim MIN MAX` | Wavelength range (nm) | _auto_ |
| `--uv_ylim YMIN YMAX` | UV-Vis y-axis limits | _auto_ |
| `--ecd_ylim YMIN YMAX` | ECD y-axis limits | _auto_ |
| `--stick` | Overlay stick spectrum | _off_ |
| `--flip_x` | Plot long-to-short wavelength | _off_ |
| `--scale FAC` | Multiply intensities after weighting | `1.0` |
| `--no_title` | Suppress figure titles | _off_ |

**Output:** `<TAG>/06_spectra/<TAG>_uvvis.png/.pdf/.csv`, `<TAG>/06_spectra/<TAG>_ecd.png/.pdf/.csv`

**Examples:**

```bash
# Single functional -- Boltzmann-weighted ECD for ephedrine
python3 or_ecd_uvvis_tools.py ephedrine/05_ecd \
    --bw ephedrine/04_boltzmann/ephedrine_energies.dat \
    --outdir ephedrine --stick --xlim 190 350

# Single functional -- UV-Vis benchmark for PNA with long-to-short axis
python3 or_ecd_uvvis_tools.py pna/05_ecd \
    --bw pna/04_boltzmann/pna_energies.dat \
    --outdir pna --stick --flip_x --xlim 200 500
```

The `**` recursive glob requires **bash >= 4.0** with `globstar` enabled (`shopt -s globstar`).
Most Linux/HPC systems meet this requirement. macOS ships bash 3.2 by default -- use
`brew install bash` or zsh (where globstar is on by default).

```bash
# Recursive glob -- collect all logs under 05_ecd/ regardless of subdirectory depth;
# use --prefix to name the output files explicitly
shopt -s globstar   # bash only; omit in zsh
python3 or_ecd_uvvis_tools.py pna/05_ecd/**/*.log \
    --bw pna/04_boltzmann/pna_energies.dat \
    --prefix pna/06_spectra/pna_wb97xd3_cpcm \
    --title "PNA / wB97X-D3 / def2-TZVP / CPCM(water)" \
    --flip_x --xlim 200 500
```

### or_vcd_ir_tools.py -- Vibrational Spectral Broadening & Plotting

```
or_vcd_ir_tools.py [OPTIONS] LOGS [LOGS ...]
```

`LOGS` can be individual `.log` files, glob patterns, or directories (searched recursively).

| Flag | Description | Default |
|------|-------------|---------|
| `--bw PATH` | Boltzmann weight file (`_energies.dat` from Stage 5) | _equal weights_ |
| `--outdir TAG` | Molecule directory; outputs go to `<TAG>/06_spectra/<TAG>_*` | |
| `--prefix STR` | Explicit output prefix (overrides `--outdir`) | `vib` |
| `--title STR` | Plot title (overrides auto-generated title) | _derived from prefix_ |
| `--ir_fwhm CM` | Gaussian FWHM for IR (cm-1) | `10` |
| `--vcd_fwhm CM` | Gaussian FWHM for VCD (cm-1) | `6` |
| `--xlim MIN MAX` | Wavenumber range (cm-1) | _auto_ |
| `--ir_ylim YMIN YMAX` | IR y-axis limits | _auto_ |
| `--vcd_ylim YMIN YMAX` | VCD y-axis limits | _auto_ |
| `--stick` | Overlay stick spectrum | _off_ |
| `--invert_ir` | Plot IR absorption peaks downward | _off_ |
| `--no_title` | Suppress figure titles | _off_ |

**Output:** `<TAG>/06_spectra/<TAG>_ir.png/.pdf/.csv`, `<TAG>/06_spectra/<TAG>_vcd.png/.pdf/.csv`

**Examples:**

```bash
# Single functional -- Boltzmann-weighted IR/VCD for methyloxirane
python3 or_vcd_ir_tools.py methyloxirane/05_vcd \
    --bw methyloxirane/04_boltzmann/methyloxirane_energies.dat \
    --outdir methyloxirane --stick --xlim 800 3200

# Single functional -- adjusted broadening for comparison with experiment
python3 or_vcd_ir_tools.py ephedrine/05_vcd \
    --bw ephedrine/04_boltzmann/ephedrine_energies.dat \
    --outdir ephedrine --vcd_fwhm 8 --ir_fwhm 12
```

The `**` recursive glob requires **bash >= 4.0** with `globstar` enabled (`shopt -s globstar`).
Most Linux/HPC systems meet this requirement. macOS ships bash 3.2 by default -- use
`brew install bash` or zsh (where globstar is on by default).

```bash
# Recursive glob -- collect all logs under 05_vcd/ regardless of subdirectory depth;
# use --prefix to name the output files explicitly
shopt -s globstar   # bash only; omit in zsh
python3 or_vcd_ir_tools.py methyloxirane/05_vcd/**/*.log \
    --bw methyloxirane/04_boltzmann/methyloxirane_energies.dat \
    --prefix methyloxirane/06_spectra/methyloxirane_b3lyp_cpcm \
    --title "methyloxirane / B3LYP / def2-TZVP / CPCM(water)" \
    --stick --xlim 800 3200
```

---

## Supported ORCA Options

### Density Functionals

The `--method` flag accepts any ORCA-recognized functional keyword. Common choices for spectroscopy:

| Functional | Type | HF Exchange | Typical Use |
|------------|------|-------------|-------------|
| B3LYP | Hybrid GGA | 20% | Geometry optimization, IR/VCD frequencies |
| PBE0 | Hybrid GGA | 25% | Slightly better than B3LYP for many properties |
| CAM-B3LYP | Range-separated | 19-65% | Charge-transfer excitations, ECD |
| wB97X-D3 | Range-separated | 22-100% | TD-DFT benchmark standard |
| M06-2X | Hybrid meta-GGA | 54% | Main-group thermochemistry |
| r2SCAN | Meta-GGA | 0% | Modern general-purpose |

**Note on VCD:** For vibrational spectra, B3LYP is the most extensively benchmarked functional and is generally recommended as a starting point. Range-separated hybrids (CAM-B3LYP, wB97X-D3) offer no systematic advantage for harmonic frequencies and are significantly more expensive for analytic Hessian calculations.

### Basis Sets

The `--basis` flag accepts one or more ORCA basis set keywords (space-separated, quoted). The auxiliary basis (e.g., `def2/J`) enables the RI-J approximation for faster Coulomb fitting.

| Basis | Quality | Typical Use |
|-------|---------|-------------|
| `def2-SVP def2/J` | Double-zeta | Geometry optimization |
| `def2-TZVP def2/J` | Triple-zeta | TD-DFT and frequency production runs |
| `def2-TZVPP def2/J` | Triple-zeta + extra polarization | Basis set convergence checks |
| `def2-QZVP def2/J` | Quadruple-zeta | Near-complete basis limit |

### Solvents

The `--solvent` flag accepts ORCA solvent keywords for the SMD implicit solvation model. Common options:

`water`, `methanol`, `ethanol`, `acetone`, `acetonitrile`, `ch2cl2`, `chcl3`, `benzene`, `toluene`, `hexane`, `cyclohexane`, `dmso`, `dmf`, `thf`, `pyridine`, `diethylether`

### Dispersion Corrections

The `--disp` flag controls how dispersion is applied:

| Value | Behavior |
|-------|-----------|
| `auto` _(default)_ | Detects whether the functional already includes dispersion (e.g., wB97X-**D3**, B97-**D**). If not, adds D3BJ. |
| `none` | No dispersion correction. |
| `D3BJ` | Grimme's D3 with Becke-Johnson damping. |
| `D4` | Grimme's D4 dispersion. |

---

## Charge and Multiplicity

All scripts parse charge and multiplicity from line 2 of the XYZ file. Two formats are supported:

```
14                           # atom count
charge=0 mult=1              # key=value (any order)
C  0.000  0.000  0.000
...
```

```
14                           # atom count
0 1                          # bare integers: charge multiplicity
C  0.000  0.000  0.000
...
```

If neither format is detected, the defaults `charge=0 mult=1` are used.

---

## Molecule List Files

Instead of passing molecule tags as positional arguments, you can supply a text file with `--list`:

```
# molecules.txt
# Lines starting with # are ignored
pna
ephedrine
aspirin
```

```bash
./1-orca-init-opt.sh --local --cpus 4 --list molecules.txt
```

---

## Example Molecules

The `pre_xyz/` directory includes starting geometries for several test molecules:

**Session 1 - Rigid chromophores (UV-Vis benchmarking):**
- `pna.xyz` - para-nitroaniline (intramolecular charge transfer)
- `dmabn.xyz` - 4-(dimethylamino)benzonitrile (dual fluorescence)
- `azobenzene.xyz` - trans-azobenzene (n->pi* and pi->pi* transitions)

**Session 2 - Flexible chiral molecules (ECD + conformational averaging):**
- `ephedrine.xyz` - (1R,2S)-ephedrine (amino alcohol, 4-10 conformers)
- `norephedrine.xyz` - (1S,2R)-norephedrine (amino alcohol)
- `methyloxirane.xyz` - (R)-methyloxirane (propylene oxide, small chiral reference)
- `phenylethan1ol.xyz` - (S)-1-phenylethan-1-ol (alpha-methylbenzyl alcohol)
- `sparteine.xyz` - (-)-sparteine (tetracyclic alkaloid, large conformer space)

**Session 3 - Vibrational chiroptical spectroscopy (VCD):**
- `methyloxirane.xyz` - (R)-methyloxirane (the canonical VCD benchmark, rigid)
- `phenylethan1ol.xyz` - (S)-1-phenylethan-1-ol (moderate flexibility, VCD literature reference)

---

## Lab Sessions

This workflow supports three independent teaching sessions. Each lab is self-contained with its own handout and report template; instructors can offer them in any order or combination.

| Lab | Handout | Report Template | Duration |
|-----|---------|-----------------|----------|
| UV-Vis Benchmarking (PNA) | `Session1_UVVis_Benchmarking_Handout` | `Student_Report_Template_UVVis_ECD` | 3 hours |
| ECD + Conformational Averaging (Ephedrine) | `Session2_ECD_Conformational_Averaging_Handout` | `Student_Report_Template_UVVis_ECD` | 3 hours |
| VCD/IR Spectroscopy (Methyloxirane) | `Session_VCD_IR_Handout` | `Student_Report_Template_VCD` | 3 hours |

The UV-Vis and ECD labs share a report template because they both use the electronic spectra branch of the pipeline (TD-DFT, then `or_ecd_uvvis_tools.py`). The VCD lab has its own dedicated template covering vibrational-specific topics (harmonic scaling, basis set convergence of VCD signs).

---

## Troubleshooting

**ORCA not found or wrong version used:**  The recommended setup is to configure `cluster.cfg` (see Quick Start step 5). This is the most reliable approach because it is explicit, per-repo, and never conflicts with system `PATH` or shell rc files. Alternatives in order of priority:
```bash
# Option 1 -- cluster.cfg (recommended, set once per clone)
echo 'ORCA_BIN=/path/to/orca' >> cluster.cfg

# Option 2 -- environment variable
export ORCA_BIN=/path/to/orca

# Option 3 -- per-invocation flag
./1-orca-init-opt.sh --orca-bin /path/to/orca pna
```
The scripts resolve the binary in this order: `--orca-bin` flag, then `ORCA_BIN` env var (including from `cluster.cfg`), then `which orca`. If the wrong version is being picked up, check each level. A common pitfall is an old ORCA installation left on `PATH` in `~/.bashrc`, which will be silently used if neither `ORCA_BIN` nor `--orca-bin` is set.

**ORCA binary path is wrong or not executable:**  The scripts validate the resolved binary before submitting any job. If you see `does not exist or is not executable`, confirm the full path in `cluster.cfg` is correct and that the binary has execute permission (`ls -l $ORCA_BIN`).

**OpenMPI version mismatch:**  ORCA 6 shared builds require the specific OpenMPI version they were compiled against to be on `LD_LIBRARY_PATH`. Set `OMPI_DIR` in `cluster.cfg` to the matching OpenMPI root; the SLURM job scripts will add `$OMPI_DIR/lib` automatically. See [INSTALL.md](INSTALL.md) for details.

**Convergence failure:** Increase `--max-iter` or try a different initial geometry. Check the `.log` file for SCF convergence messages.

**Empty Boltzmann output:** All conformers fell below `--p-cut`. Lower the threshold (e.g., `--p-cut 0.001`) or check that Stage 4 completed successfully for all conformers.

**Plotting tool import errors:** Install Python dependencies: `pip install -r requirements.txt`

**CREST fails with a shared library error on HPC** (e.g., `version 'GLIBCXX_3.4.XX' not found`): The CREST binary was compiled against a newer libgcc/libgomp than the cluster's system provides. Create a minimal conda environment to supply the missing libraries, then activate it before submitting the Stage 2b job:

```bash
conda create -n xtb_rt -c conda-forge libgcc libgomp libgfortran
conda activate xtb_rt
bash 2b-crest-conf-search.sh TAG
```

The SLURM job inherits the activated environment's LD_LIBRARY_PATH, so no changes to the batch script are needed. See [INSTALL.md](INSTALL.md) Section 6 for full setup details.

**VCD frequency job runs slowly:** Analytic Hessian calculations scale more steeply than single-point energies. Use `--dry-run` first to check conformer count, and reduce cores only if memory is the bottleneck. B3LYP is recommended over range-separated hybrids for frequency calculations.

---

## Citation

If you use this workflow in published research or teaching, please cite:

> [Author], [Author]. A Computational Spectroscopy Module for Teaching DFT-Based Spectral Workflows in Graduate Chemistry. _J. Chem. Educ._ **2026**, _XX_, XXXX-XXXX.

---

## License

[To be determined - recommend MIT or BSD-3-Clause for JCE educational software.]
