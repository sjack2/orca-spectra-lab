# Installation Guide

This guide covers installing all software required to run the ORCA Electronic Spectra Lab workflow.  Three environments are supported:

- **Linux workstation** (native Ubuntu, Fedora, or similar)
- **Windows via WSL** (Windows Subsystem for Linux)
- **HPC cluster** (SLURM-managed, modules-based)

---

## Prerequisites at a Glance

| Software | Version | Required? | Purpose |
|----------|---------|-----------|---------|
| ORCA | ≥ 5.0 | Yes | Quantum chemistry calculations |
| OpenMPI | ≥ 4.0 | Yes | Parallel execution for ORCA |
| Open Babel | ≥ 3.0 | Yes | XYZ/SDF conversion, Confab conformer search |
| Python 3 | ≥ 3.8 | Yes | Spectral plotting tool |
| NumPy | ≥ 1.20 | Yes | Numerical operations |
| Matplotlib | ≥ 3.4 | Yes | Plot generation |
| Pandas | ≥ 1.3 | Yes | Data handling |
| CREST | ≥ 2.12 | Optional | GFN2-xTB conformer search (Stage 2b) |
| Bash | ≥ 4.0 | Yes | Workflow scripts (ships with all Linux distributions) |
| getopt (GNU) | — | Yes | CLI flag parsing (ships with `util-linux`) |

---

## 1. Check Existing Installations

Before installing anything, check what you already have:

```bash
which orca        # ORCA quantum chemistry
orca --version    # (some builds print version info)
which obabel      # Open Babel
obabel -V         # Open Babel version
python3 --version # Python 3
mpirun --version  # OpenMPI
which crest       # CREST (optional)
```

If a command prints a valid path and version, that tool is already installed and you can skip its section below.

---

## 2. Python 3

Most Linux distributions ship Python 3.  Verify:

```bash
python3 --version
```

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
```

### Fedora / CentOS / RHEL

```bash
sudo dnf install python3 python3-pip
```

### Python Dependencies

Install the plotting tool's requirements:

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install numpy matplotlib pandas
```

> **Tip:** If you prefer isolated environments, use a virtual environment:
> ```bash
> python3 -m venv .venv
> source .venv/bin/activate
> pip install -r requirements.txt
> ```

---

## 3. Open Babel

Open Babel provides the `obabel` command for file format conversion and the Confab conformer search engine used in Stage 2.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install openbabel
```

### Fedora / CentOS

```bash
sudo dnf install openbabel
```

### Verify

```bash
obabel -V
```

This should print the version (e.g., `Open Babel 3.1.1`).  Test Confab:

```bash
obabel -L conformer
```

This should list the Confab conformer search plugin.

---

## 4. OpenMPI

ORCA uses OpenMPI for parallel (multi-core) calculations.  **Important:** the OpenMPI version must be compatible with the ORCA build you download.  ORCA 6.0 ships with OpenMPI 4.1.x; ORCA 5.0 typically expects OpenMPI 4.1.x as well.  Check the ORCA release notes for the exact required version.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install openmpi-bin libopenmpi-dev
```

### Fedora / CentOS

```bash
sudo dnf install openmpi openmpi-devel
```

On Fedora/CentOS, you may need to activate the MPI environment:

```bash
module load mpi/openmpi-x86_64   # or source /etc/profile.d/modules.sh first
```

### Verify

```bash
mpirun --version
```

Should print something like `mpirun (Open MPI) 4.1.4`.

### Version Mismatch Troubleshooting

If ORCA produces errors like `mca_base_component_find: unable to open ...` or segfaults during parallel runs, the OpenMPI version does not match what ORCA expects.  Solutions:

1. **Use ORCA's bundled MPI** — some ORCA distributions include OpenMPI libraries.  Add them to `LD_LIBRARY_PATH`:
   ```bash
   export LD_LIBRARY_PATH=/opt/orca_6_0_1:$LD_LIBRARY_PATH
   ```

2. **Install the matching version from source** — download from [open-mpi.org](https://www.open-mpi.org/software/) and build:
   ```bash
   tar -xf openmpi-4.1.6.tar.gz
   cd openmpi-4.1.6
   ./configure --prefix=$HOME/openmpi
   make -j$(nproc) && make install
   export PATH=$HOME/openmpi/bin:$PATH
   export LD_LIBRARY_PATH=$HOME/openmpi/lib:$LD_LIBRARY_PATH
   ```

---

## 5. ORCA

ORCA is free for academic use but requires registration.

### a) Register and Download

1. Go to the [ORCA Forum](https://orcaforum.kofo.mpg.de).
2. Register for an account.
3. Navigate to the Downloads section and download the latest Linux version (e.g., `orca_6_0_1_linux_x86-64_shared_openmpi416.tar.xz`).
4. Choose the **shared library** build that matches your OpenMPI version.

### b) Extract and Install

```bash
# Extract
tar -xf orca_6_0_1_linux_x86-64_shared_openmpi416.tar.xz

# Move to a permanent location
sudo mv orca_6_0_1 /opt/orca_6_0_1

# Or, for a user-local install (no sudo required):
mv orca_6_0_1 $HOME/orca_6_0_1
```

### c) Configure PATH

Add to your `~/.bashrc` (or `~/.bash_profile`, `~/.zshrc`, etc.):

```bash
# ORCA
export PATH="/opt/orca_6_0_1:$PATH"
export LD_LIBRARY_PATH="/opt/orca_6_0_1:$LD_LIBRARY_PATH"
```

Then reload:

```bash
source ~/.bashrc
```

Alternatively, the workflow scripts accept `--orca-bin /path/to/orca` or respect the `ORCA_BIN` environment variable:

```bash
export ORCA_BIN=/opt/orca_6_0_1/orca
```

### d) Verify

```bash
which orca

# Run a quick test (single-point HF on water)
mkdir -p /tmp/orca_test && cd /tmp/orca_test
cat > test.inp << 'EOF'
! HF STO-3G
* xyz 0 1
O   0.0   0.0   0.0
H   0.0   0.757  0.587
H   0.0  -0.757  0.587
*
EOF
orca test.inp > test.out 2>&1
grep "FINAL SINGLE POINT ENERGY" test.out
# Should print: FINAL SINGLE POINT ENERGY      -74.96...
cd - && rm -rf /tmp/orca_test
```

---

## 6. CREST (Optional)

CREST performs conformer searches using GFN2-xTB metadynamics.  It is used in the alternative Stage 2b path and is more thorough than Confab for flexible molecules.

### Download

CREST binaries are available from the [Grimme group GitHub](https://github.com/crest-lab/crest/releases).

```bash
# Download the latest static binary
wget https://github.com/crest-lab/crest/releases/download/v3.0.2/crest-gnu-12-ubuntu-latest.tar.xz
tar -xf crest-gnu-12-ubuntu-latest.tar.xz
sudo mv crest /usr/local/bin/
```

Or install via conda:

```bash
conda install -c conda-forge crest
```

### Verify

```bash
crest --version
```

### Runtime libraries on HPC

On some HPC clusters the system's default libgcc and libgomp are older than
what the CREST binary was compiled against.  The symptom is an error like
`version 'GLIBCXX_3.4.XX' not found` when you run `crest --version`.

The cleanest fix is a minimal conda environment that supplies up-to-date
runtime libraries without installing CREST itself through conda:

```bash
conda create -n xtb_rt -c conda-forge libgcc libgomp libgfortran
conda activate xtb_rt
```

Place the CREST binary somewhere on your PATH, for example `~/opt/bin/`:

```bash
mkdir -p ~/opt/bin
cp crest ~/opt/bin/crest
chmod +x ~/opt/bin/crest
echo 'export PATH=$HOME/opt/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Activate the environment before submitting any SLURM job that calls CREST:

```bash
conda activate xtb_rt
bash 2b-crest-conf-search.sh TAG
```

The SLURM job inherits the submitting shell's PATH and LD_LIBRARY_PATH, so
`crest` remains findable on the compute node without any extra module
commands in the batch script.

---

## 7. Windows Users: Installing WSL

If you are on Windows 10 (version 2004+) or Windows 11, you can run the entire workflow inside the Windows Subsystem for Linux.

### a) Check Windows Version

Press `Win+R`, type `winver`, and press Enter.  Confirm you have Windows 10 version 2004 or later, or Windows 11.

### b) Install WSL

Open **PowerShell** or **Command Prompt** as Administrator and run:

```powershell
wsl --install
```

This enables WSL, installs the Virtual Machine Platform, and installs Ubuntu by default.  Restart your computer when prompted.

If the automatic method does not work, enable the features manually:

```powershell
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```

Then restart, and install Ubuntu from the Microsoft Store.

### c) First Launch

Launch Ubuntu from the Start Menu.  You will be prompted to create a Linux username and password.

### d) Verify

In PowerShell:

```powershell
wsl --list --verbose
```

This should show your Ubuntu installation with State = Running.

### e) Install Software Inside WSL

Once inside the WSL terminal, follow Sections 2–6 above as if you were on a native Linux machine.  All `apt` and `pip` commands work identically.

> **Performance note:** WSL2 provides near-native performance for CPU-bound tasks like ORCA.  Store your working files inside the Linux filesystem (`/home/username/`) rather than on the Windows drive (`/mnt/c/`) for best I/O performance.

---

## 8. HPC Cluster Setup

On a SLURM-managed cluster, software is typically loaded via the `module` system rather than installed locally.

### Typical Module Setup

```bash
module load apps/orca/6.0.1
module load compilers/gcc/12.2.0
module load mpi/openmpi/4.1.6
module load apps/python/3.11
```

The exact module names vary by institution.  Check with `module avail orca` and `module avail openmpi`.

### Python Dependencies

On clusters, you may not have root access.  Use `--user` or a virtual environment:

```bash
pip install --user -r requirements.txt
# or
python3 -m venv ~/.venvs/spectra
source ~/.venvs/spectra/bin/activate
pip install -r requirements.txt
```

### SLURM Configuration

The workflow scripts accept SLURM options:

```bash
./1-orca-init-opt.sh --partition main --time 02:00:00 --cpus 8 --mem-per-cpu 4096 aspirin
```

The default partition is `circe` (University of South Florida's cluster).  Override with `--partition` to match your institution.

### CREST on HPC

CREST is often available as a module:

```bash
module load apps/crest/3.0
```

If not, download the static binary and place it in your `$HOME/bin`
(or `$HOME/opt/bin/` -- any directory on your PATH):

```bash
mkdir -p $HOME/opt/bin
wget -O $HOME/opt/bin/crest https://github.com/crest-lab/crest/releases/download/v3.0.2/crest-gnu-12-ubuntu-latest.tar.xz
tar -xf crest-gnu-12-ubuntu-latest.tar.xz
mv crest $HOME/opt/bin/crest
chmod +x $HOME/opt/bin/crest
export PATH=$HOME/opt/bin:$PATH
```

If the cluster's system libraries are too old for the binary (a
`GLIBCXX version not found` error), see the "Runtime libraries on HPC"
subsection in Section 6 above for the conda-based fix.

---

## 9. Final Verification

Run these commands to confirm everything is ready:

```bash
echo "=== Python ==="
python3 --version
python3 -c "import numpy, matplotlib, pandas; print('NumPy', numpy.__version__); print('Matplotlib', matplotlib.__version__); print('Pandas', pandas.__version__)"

echo "=== Open Babel ==="
obabel -V

echo "=== OpenMPI ==="
mpirun --version | head -1

echo "=== ORCA ==="
which orca
echo "(run 'orca --version' or a test job to confirm)"

echo "=== CREST (optional) ==="
which crest 2>/dev/null && crest --version || echo "CREST not installed (optional)"

echo "=== Workflow scripts ==="
./1-orca-init-opt.sh --help | head -3
```

If all commands succeed, you are ready to run the pipeline.  See [README.md](README.md) for usage instructions.
