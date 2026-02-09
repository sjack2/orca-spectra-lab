# ORCA Electronic Spectra Lab

A modular Bash/Python workflow for computing UV-Vis absorption and electronic circular dichroism (ECD) spectra using ORCA time-dependent density functional theory (TD-DFT).  Designed for teaching and research — runs identically on a local workstation or a SLURM-managed HPC cluster.

> **Accompanying article:**  _A Computational Spectroscopy Module for Teaching TD-DFT Workflows in Graduate Chemistry_ (in preparation for the _Journal of Chemical Education_).

---

## Features

- **Six-stage pipeline** from raw XYZ geometry to publication-quality broadened spectra.
- **Dual execution mode** — every ORCA-calling script auto-detects whether SLURM is available.  Pass `--local` to force direct execution; omit it on a cluster and jobs are submitted via `sbatch`.
- **Two conformer search paths** — Open Babel Confab (fast, force-field-based) or CREST (GFN2-xTB metadynamics, more thorough).
- **Boltzmann-weighted spectral averaging** with configurable temperature and population threshold.
- **Physically correct broadening** — Gaussian convolution in energy space (eV) with Jacobian correction for ECD, producing both PNG/PDF plots and CSV data tables.
- **`--dry-run` on every script** — inspect generated ORCA inputs without running any calculations.

---

## Pipeline Overview

```
Stage 1                    Gas-phase geometry optimisation
  │                        1-orca-init-opt.sh
  ▼
Stage 2                    Conformer enumeration
  │                        2-orca-conf-search.sh  (Confab)
  │                   or   2b-crest-conf-search.sh (CREST)
  ▼
Stage 3                    Split conformers into individual XYZ files
  │                        3-orca-conf-split.sh   (Confab output)
  │                   or   3b-crest-conf-split.sh  (CREST output)
  ▼
Stage 4                    Solvent-phase re-optimisation (SMD/CPCM)
  │                        4-orca-solvent-opt.sh
  ▼
Stage 5                    Boltzmann weighting & filtering
  │                        5-orca-boltzmann-weight.sh
  ▼
Stage 6                    TD-DFT excited-state calculations
  │                        6-orca-ecd.sh
  ▼
Plot                       Broadened UV-Vis & ECD spectra
                           or_ecd_uvvis_tools.py
```

---

## Quick Start

### 1. Install prerequisites

See [INSTALL.md](INSTALL.md) for detailed instructions.  In brief, you need:

| Software | Version | Purpose |
|----------|---------|---------|
| ORCA | ≥ 5.0 | Quantum chemistry engine |
| OpenMPI | ≥ 4.0 | Parallel execution for ORCA |
| Open Babel | ≥ 3.0 | File conversion & Confab conformer search |
| Python 3 | ≥ 3.8 | Plotting tool |
| CREST | ≥ 2.12 | _(optional)_ GFN2-xTB conformer search |

### 2. Clone the repository

```bash
git clone https://github.com/sjack2/orca-electronic-spectra-lab.git
cd orca-electronic-spectra-lab
```

### 3. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 4. Make scripts executable

```bash
chmod +x *.sh
```

### 5. Set up your molecule

Place a starting XYZ file in `pre_xyz/`.  Line 2 must encode the charge and multiplicity:

```
14
charge=0 mult=1
C    0.000000    0.000000    0.000000
...
```

Both `charge=0 mult=1` (key=value) and `0 1` (bare integers) are accepted.  If neither is found, defaults are charge=0, mult=1.

### 6. Run the pipeline (local workstation example)

```bash
# Stage 1: Optimise geometry in vacuum
./1-orca-init-opt.sh --local --cpus 4 pna

# Stage 2+3: Generate and split conformers (Confab path)
./2-orca-conf-search.sh --ecut 5 --conf 500 pna
./3-orca-conf-split.sh pna

# Stage 4: Re-optimise each conformer in solvent
./4-orca-solvent-opt.sh --local --cpus 4 --solvent water pna

# Stage 5: Boltzmann filter
./5-orca-boltzmann-weight.sh pna

# Stage 6: TD-DFT on populated conformers
./6-orca-ecd.sh --local --cpus 4 --method wB97X-D3 --roots 20 pna

# Plot UV-Vis spectrum
python3 or_ecd_uvvis_tools.py pna/ecd \
    --bw pna/bw_results/pna_energies.dat \
    --prefix pna_uv --xlim 200 500
```

### 7. Run the pipeline (HPC/SLURM example)

On a cluster with SLURM, simply omit `--local`.  The scripts detect `sbatch` automatically and submit jobs:

```bash
./1-orca-init-opt.sh --cpus 8 --partition main --time 02:00:00 pna
# → submits SLURM job; monitor with squeue
```

---

## Directory Structure

Every script follows the same directory convention.  For a molecule tagged `aspirin`:

```
aspirin/
├── aspirin_orca_opt/          Stage 1  — gas-phase optimisation
│   ├── aspirin.inp
│   ├── aspirin.log
│   └── aspirin.xyz            optimised geometry
├── orca_opt_conf/             Stage 2+3 — conformer search
│   ├── aspirin.sdf
│   ├── aspirin_combined.sdf
│   ├── split_sdf/
│   └── split_xyz/             individual conformer XYZ files
│       ├── aspirin_1.xyz
│       ├── aspirin_2.xyz
│       └── ...
├── solvent_opt/               Stage 4  — solvent-phase optimisation
│   ├── aspirin_1/
│   │   ├── aspirin_1.inp
│   │   ├── aspirin_1.log
│   │   └── aspirin_1.xyz
│   └── aspirin_2/
│       └── ...
├── bw_results/                Stage 5  — Boltzmann weighting
│   ├── aspirin_energies.dat   full table (CID, E, ΔE, p)
│   └── aspirin_bw_labels.dat  conformer IDs above threshold
└── ecd/                       Stage 6  — TD-DFT
    ├── aspirin_1/
    │   ├── aspirin_1.inp
    │   └── aspirin_1.log
    └── aspirin_2/
        └── ...
```

Starting geometries live in `pre_xyz/`:

```
pre_xyz/
├── aspirin.xyz
├── pna.xyz
└── ephedrine.xyz
```

---

## CLI Reference

All scripts accept `--help` for full usage.  Flags shown with `[default]`.

### 1-orca-init-opt.sh — Gas-Phase Geometry Optimisation

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
| `--grid N` | `-g` | DEFGRID level (1–3) | `3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--orca-bin PATH` | | Path to ORCA binary | _auto-detected_ |
| `--list FILE` | | Text file of TAGs (one per line) | |
| `--local` | | Run ORCA directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `circe` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Output:** `<TAG>/<TAG>_orca_opt/<TAG>.xyz` (optimised geometry)

### 2-orca-conf-search.sh — Confab Conformer Search

```
2-orca-conf-search.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--conf N` | Maximum conformers | `1000` |
| `--ecut E` | Energy cutoff (kcal/mol) | `6` |
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo commands without running | |

**Output:** `<TAG>/orca_opt_conf/<TAG>_combined.sdf`

### 2b-crest-conf-search.sh — CREST Conformer Search

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
| `--partition NAME` | | SLURM partition _(HPC only)_ | `circe` |
| `--time HH:MM:SS` | | SLURM wall time _(HPC only)_ | `06:00:00` |

**Output:** `<TAG>/orca_opt_conf/crest_conformers.xyz`

### 3-orca-conf-split.sh — Split Confab Output

```
3-orca-conf-split.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo actions without running | |

**Input:** `<TAG>/orca_opt_conf/<TAG>_combined.sdf`
**Output:** `<TAG>/orca_opt_conf/split_xyz/<TAG>_N.xyz`

### 3b-crest-conf-split.sh — Split CREST Output

```
3b-crest-conf-split.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo actions without running | |

**Input:** `<TAG>/orca_opt_conf/crest_conformers.xyz`
**Output:** `<TAG>/orca_opt_conf/split_xyz/<TAG>_NNN.xyz` (zero-padded)

### 4-orca-solvent-opt.sh — Solvent-Phase Re-optimisation

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
| `--grid N` | `-g` | DEFGRID level (1–3) | `3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--orca-bin PATH` | | Path to ORCA binary | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run ORCA directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `circe` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Input:** `<TAG>/orca_opt_conf/split_xyz/` (from Stage 3/3b)
**Output:** `<TAG>/solvent_opt/<CID>/<CID>.xyz` (optimised geometry per conformer)

### 5-orca-boltzmann-weight.sh — Boltzmann Weighting & Filtering

```
5-orca-boltzmann-weight.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--temp K` | Temperature in Kelvin | `298.15` |
| `--p-cut VAL` | Probability cutoff (e.g., 0.01 = 1%) | `0.01` |
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Show what would be computed | |

**Input:** `<TAG>/solvent_opt/<CID>/<CID>.log` (ORCA outputs from Stage 4)
**Output:**
- `<TAG>/bw_results/<TAG>_energies.dat` — full table: conformer ID, energy (Hartree), relative energy (kcal/mol), Boltzmann probability
- `<TAG>/bw_results/<TAG>_bw_labels.dat` — conformer IDs above the probability cutoff

### 6-orca-ecd.sh — TD-DFT Excited-State Calculations

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
| `--grid N` | `-g` | DEFGRID level (1–3) | `3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--orca-bin PATH` | | Path to ORCA binary | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run ORCA directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `circe` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Input:** `<TAG>/bw_results/<TAG>_bw_labels.dat` (from Stage 5) + `<TAG>/solvent_opt/<CID>/<CID>.xyz` (from Stage 4)
**Output:** `<TAG>/ecd/<CID>/<CID>.log` (ORCA TD-DFT output with absorption and CD spectrum blocks)

### or_ecd_uvvis_tools.py — Spectral Broadening & Plotting

```
or_ecd_uvvis_tools.py [OPTIONS] LOGS [LOGS ...]
```

`LOGS` can be individual `.log` files, glob patterns, or directories (searched recursively).

| Flag | Description | Default |
|------|-------------|---------|
| `--bw PATH` | Boltzmann weight file (`_energies.dat` from Stage 5) | _equal weights_ |
| `--prefix STR` | Output filename prefix | `spectra` |
| `--uv_fwhm EV` | Gaussian FWHM for UV-Vis (eV) | `0.35` |
| `--ecd_fwhm EV` | Gaussian FWHM for ECD (eV) | `0.25` |
| `--xlim MIN MAX` | Wavelength range (nm) | _auto_ |
| `--uv_ylim YMIN YMAX` | UV-Vis y-axis limits | _auto_ |
| `--ecd_ylim YMIN YMAX` | ECD y-axis limits | _auto_ |
| `--stick` | Overlay stick spectrum | _off_ |
| `--flip_x` | Plot long→short wavelength | _off_ |
| `--scale FAC` | Multiply intensities after weighting | `1.0` |
| `--no_title` | Suppress figure titles | _off_ |

**Output:** `<prefix>_uv.png`, `<prefix>_uv.pdf`, `<prefix>_ecd.png`, `<prefix>_ecd.pdf`, `<prefix>_uv.csv`, `<prefix>_ecd.csv`

**Examples:**

```bash
# Boltzmann-weighted UV-Vis for ephedrine
python3 or_ecd_uvvis_tools.py ephedrine/ecd \
    --bw ephedrine/bw_results/ephedrine_energies.dat \
    --prefix ephedrine_spectra --xlim 190 320

# Single-molecule UV-Vis with stick overlay
python3 or_ecd_uvvis_tools.py pna/ecd \
    --bw pna/bw_results/pna_energies.dat \
    --prefix pna_uv --stick --xlim 200 500

# Compare two functionals (separate runs, overlaid manually)
python3 or_ecd_uvvis_tools.py pna/ecd_B3LYP \
    --prefix pna_B3LYP --xlim 200 500
python3 or_ecd_uvvis_tools.py pna/ecd_wB97XD3 \
    --prefix pna_wB97XD3 --xlim 200 500
```

---

## Supported ORCA Options

### Density Functionals

The `--method` flag accepts any ORCA-recognised functional keyword.  Common choices for spectroscopy:

| Functional | Type | HF Exchange | Typical Use |
|------------|------|-------------|-------------|
| B3LYP | Hybrid GGA | 20% | Geometry optimisation, general-purpose |
| PBE0 | Hybrid GGA | 25% | Slightly better than B3LYP for many properties |
| CAM-B3LYP | Range-separated | 19–65% | Charge-transfer excitations, ECD |
| ωB97X-D3 | Range-separated | 22–100% | TD-DFT benchmark standard |
| M06-2X | Hybrid meta-GGA | 54% | Main-group thermochemistry |
| r2SCAN | Meta-GGA | 0% | Modern general-purpose |

### Basis Sets

The `--basis` flag accepts one or more ORCA basis set keywords (space-separated, quoted).  The auxiliary basis (e.g., `def2/J`) enables the RI-J approximation for faster Coulomb fitting.

| Basis | Quality | Typical Use |
|-------|---------|-------------|
| `def2-SVP def2/J` | Double-ζ | Geometry optimisation |
| `def2-TZVP def2/J` | Triple-ζ | TD-DFT production runs |
| `def2-TZVPP def2/J` | Triple-ζ + extra polarisation | Basis set convergence checks |
| `def2-QZVP def2/J` | Quadruple-ζ | Near-complete basis limit |

### Solvents

The `--solvent` flag accepts ORCA solvent keywords for the SMD implicit solvation model.  Common options:

`water`, `methanol`, `ethanol`, `acetone`, `acetonitrile`, `ch2cl2`, `chcl3`, `benzene`, `toluene`, `hexane`, `cyclohexane`, `dmso`, `dmf`, `thf`, `pyridine`, `diethylether`

### Dispersion Corrections

The `--disp` flag controls how dispersion is applied:

| Value | Behaviour |
|-------|-----------|
| `auto` _(default)_ | Detects whether the functional already includes dispersion (e.g., ωB97X-**D3**, B97-**D**).  If not, adds D3BJ. |
| `none` | No dispersion correction. |
| `D3BJ` | Grimme's D3 with Becke–Johnson damping. |
| `D4` | Grimme's D4 dispersion. |

---

## Charge and Multiplicity

All scripts parse charge and multiplicity from line 2 of the XYZ file.  Two formats are supported:

```
14                           ← atom count
charge=0 mult=1              ← key=value (any order)
C  0.000  0.000  0.000
...
```

```
14                           ← atom count
0 1                          ← bare integers: charge multiplicity
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

**Session 1 — Rigid chromophores (UV-Vis benchmarking):**
- `pna.xyz` — para-nitroaniline (intramolecular charge transfer)
- `dmabn.xyz` — 4-(dimethylamino)benzonitrile (dual fluorescence)
- `azobenzene.xyz` — trans-azobenzene (n→π* and π→π* transitions)

**Session 2 — Flexible chiral molecules (ECD + conformational averaging):**
- `ephedrine.xyz` — (1R,2S)-ephedrine (amino alcohol, 4–10 conformers)
- `norephedrine.xyz` — (1S,2R)-norephedrine (amino alcohol)
- `methyloxirane.xyz` — (R)-methyloxirane (propylene oxide, small chiral reference)
- `phenylethan1ol.xyz` — (S)-1-phenylethan-1-ol (α-methylbenzyl alcohol)
- `sparteine.xyz` — (−)-sparteine (tetracyclic alkaloid, large conformer space)

---

## Troubleshooting

**ORCA not found:** Set the `ORCA_BIN` environment variable or use the `--orca-bin` flag:
```bash
export ORCA_BIN=/opt/orca_6_0_1/orca
# or
./1-orca-init-opt.sh --orca-bin /opt/orca_6_0_1/orca --local pna
```

**OpenMPI version mismatch:** ORCA ships with a specific OpenMPI version.  If you get MPI errors, ensure the OpenMPI in your `PATH` matches the version ORCA was compiled against.  See [INSTALL.md](INSTALL.md) for details.

**Convergence failure:** Increase `--max-iter` or try a different initial geometry.  Check the `.log` file for SCF convergence messages.

**Empty Boltzmann output:** All conformers fell below `--p-cut`.  Lower the threshold (e.g., `--p-cut 0.001`) or check that Stage 4 completed successfully for all conformers.

**Plotting tool import errors:** Install Python dependencies: `pip install -r requirements.txt`

---

## Citation

If you use this workflow in published research or teaching, please cite:

> [Author], [Author]. A Computational Spectroscopy Module for Teaching TD-DFT Workflows in Graduate Chemistry. _J. Chem. Educ._ **2026**, _XX_, XXXX–XXXX.

---

## License

[To be determined — recommend MIT or BSD-3-Clause for JCE educational software.]
