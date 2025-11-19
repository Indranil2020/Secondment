#!/usr/bin/env python
"""Test if B3LYP is GPU-compatible"""

from pyscf import gto
from gpu4pyscf import dft as gpu_dft

# Create simple test molecule
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', verbose=0)

print("Testing B3LYP GPU support...")
try:
    mf = gpu_dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.verbose = 0
    mf.kernel()
    print("✓ B3LYP IS GPU-COMPATIBLE!")
    print(f"  Energy: {mf.e_tot:.6f} a.u.")
except Exception as e:
    print(f"✗ B3LYP NOT GPU-compatible")
    print(f"  Error: {str(e)[:200]}")

print("\nTesting PBE GPU support...")
try:
    mf = gpu_dft.RKS(mol)
    mf.xc = 'pbe'
    mf.verbose = 0
    mf.kernel()
    print("✓ PBE IS GPU-COMPATIBLE!")
    print(f"  Energy: {mf.e_tot:.6f} a.u.")
except Exception as e:
    print(f"✗ PBE NOT GPU-compatible")
    print(f"  Error: {str(e)[:200]}")
