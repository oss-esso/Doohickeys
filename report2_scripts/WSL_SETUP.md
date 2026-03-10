# WSL Setup Guide for VQE PES Explorer

## Why WSL?

PySCF (Python Simulations of Chemistry Framework) is not natively available on Windows due to its dependency on compiled Fortran/C libraries. WSL (Windows Subsystem for Linux) provides a full Linux environment that supports PySCF.

## Prerequisites

- Windows 10 version 2004+ or Windows 11
- Administrator access to install WSL

---

## Step 1: Install WSL

Open PowerShell as Administrator and run:

```powershell
wsl --install -d Ubuntu-22.04
```

This installs WSL 2 with Ubuntu 22.04. Restart your computer when prompted.

After restart, Ubuntu will open and ask you to create a username and password.

---

## Step 2: Update Ubuntu

```bash
sudo apt update && sudo apt upgrade -y
```

---

## Step 3: Install Python and Build Dependencies

```bash
# Install Python 3.11 and development tools
sudo apt install -y python3.11 python3.11-venv python3.11-dev python3-pip

# Install build dependencies for PySCF
sudo apt install -y build-essential cmake gfortran libopenblas-dev liblapack-dev

# Install OpenGL dependencies for visualization (optional)
sudo apt install -y libgl1-mesa-dev libglfw3-dev libglew-dev
```

---

## Step 4: Set Up Project Environment

```bash
# Navigate to your project (accessible via /mnt/c/...)
cd /mnt/c/Users/Edoardo/Documents/EPFL/Doohickeys

# Create virtual environment
python3.11 -m venv .venv
source .venv/bin/activate

# Upgrade pip
pip install --upgrade pip setuptools wheel
```

---

## Step 5: Install Python Dependencies

Create a `requirements.txt` file or run:

```bash
pip install numpy scipy
pip install pyscf
pip install pennylane pennylane-qiskit
pip install qutip
pip install scikit-image
pip install matplotlib
```

For visualization (requires X server):
```bash
pip install glfw PyOpenGL PyOpenGL_accelerate
```

---

## Step 6: Verify Installation

```bash
python -c "
import pyscf
from pyscf import gto, scf

# Quick H2 test
mol = gto.Mole()
mol.atom = 'H 0 0 0; H 0 0 0.74'
mol.basis = 'sto-3g'
mol.build()

mf = scf.RHF(mol)
e_hf = mf.kernel()
print(f'H2 (STO-3G) HF Energy: {e_hf:.6f} Ha')
print('PySCF installation verified!')
"
```

Expected output:
```
H2 (STO-3G) HF Energy: -1.117501 Ha
PySCF installation verified!
```

---

## Step 7: Run the VQE PES Explorer

```bash
# Activate environment
source .venv/bin/activate

# Run basic H2 scan (no visualization)
python -m report2_scripts.main --molecule h2 --bond-range "0.5,2.5,15"

# Run with benchmarking
python -m report2_scripts.main --molecule h2 --benchmark

# Run LiH scan
python -m report2_scripts.main --molecule lih --bond-range "1.0,3.0,15" --benchmark
```

---

## Optional: GUI Support (X Server)

To use the 3D visualization (`--render`), you need an X server on Windows:

1. **Install VcXsrv** (or WSLg if on Windows 11):
   - Download from: https://sourceforge.net/projects/vcxsrv/
   - Run XLaunch with: Multiple windows, Start no client, Disable access control

2. **Configure WSL** (add to `~/.bashrc`):
   ```bash
   export DISPLAY=$(cat /etc/resolv.conf | grep nameserver | awk '{print $2}'):0.0
   export LIBGL_ALWAYS_INDIRECT=1
   ```

3. **Test**:
   ```bash
   source ~/.bashrc
   glxinfo | grep "OpenGL"
   ```

4. **Run with visualization**:
   ```bash
   python -m report2_scripts.main --molecule h2 --render --orbital 0
   ```

---

## Troubleshooting

### PySCF import error
```bash
pip uninstall pyscf
pip install pyscf --no-cache-dir
```

### OpenGL errors
```bash
export MESA_GL_VERSION_OVERRIDE=3.3
```

### Memory issues on large molecules
Increase WSL memory in `%UserProfile%\.wslconfig`:
```
[wsl2]
memory=8GB
swap=4GB
```

---

## Example Session

```bash
# Start WSL
wsl

# Navigate to project
cd /mnt/c/Users/Edoardo/Documents/EPFL/Doohickeys

# Activate environment
source .venv/bin/activate

# Run full benchmark on H2
python -m report2_scripts.main \
    --molecule h2 \
    --basis sto-3g \
    --bond-range "0.4,3.0,20" \
    --ansatz uccsd \
    --optimizer Adam \
    --maxiter 200 \
    --benchmark \
    --output output/h2_benchmark

# Results saved to output/h2_benchmark/pes_h2_sto-3g.json
```

---

## Quick Reference

| Command | Description |
|---------|-------------|
| `wsl` | Start Ubuntu |
| `source .venv/bin/activate` | Activate Python environment |
| `python -m report2_scripts.main --help` | Show all options |
| `--molecule h2/lih/h2o` | Select molecule |
| `--benchmark` | Enable multi-method comparison |
| `--render` | Enable 3D visualization |
| `exit` | Return to Windows |
