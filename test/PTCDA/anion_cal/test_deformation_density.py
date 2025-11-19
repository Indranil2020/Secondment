#!/usr/bin/env python3
"""
Test script for calculating deformation density (density difference between 
SCF molecular density and superposition of atomic densities).

This calculates:
1. Ground state SCF density (with electron-electron interaction)
2. Promolecule density (superposition of non-interacting atomic densities)
3. Deformation density = SCF density - Promolecule density

Tests both RKS and UKS for CPU and GPU.

NOTE: For GPU tests to work, you need CUDA libraries in LD_LIBRARY_PATH.
Run with: export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH
Or simply: ./run_cal.sh (which sets this automatically)
"""

import os
import sys

# Set CUDA library path for GPU4PySCF (must be done before importing)
cuda_path = '/usr/local/cuda-12.9/lib64'
if os.path.exists(cuda_path):
    current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
    if cuda_path not in current_ld_path:
        os.environ['LD_LIBRARY_PATH'] = f"{cuda_path}:{current_ld_path}"
        print(f"✓ Added CUDA path to LD_LIBRARY_PATH: {cuda_path}")

import numpy as np
from pyscf import gto, dft, scf
from pyscf.tools import cubegen

# Try to import GPU4PySCF
try:
    from gpu4pyscf import dft as gpu_dft
    GPU_AVAILABLE = True
    print("✓ GPU4PySCF available")
except ImportError as e:
    GPU_AVAILABLE = False
    print(f"✗ GPU4PySCF not available: {e}")
    print("  Hint: Make sure CUDA libraries are in LD_LIBRARY_PATH")
    print("  Run: export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH")
except Exception as e:
    GPU_AVAILABLE = False
    print(f"✗ GPU4PySCF import error: {type(e).__name__}: {e}")
    print("  This might be a CUDA runtime issue")

def test_deformation_density_rks_cpu():
    """Test deformation density for RKS (closed-shell) on CPU"""
    print("\n" + "="*70)
    print("TEST 1: RKS (Closed-shell) - CPU")
    print("="*70)
    
    # Create H2O molecule (neutral, closed-shell)
    mol = gto.M(
        atom='O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692',
        basis='6-31g',
        charge=0,
        spin=0,
        verbose=0
    )
    
    print(f"Molecule: H2O (neutral)")
    print(f"Electrons: {mol.nelectron}")
    print(f"Basis functions: {mol.nao_nr()}")
    
    # 1. Calculate SCF density (with electron-electron interaction)
    print("\n1. Calculating SCF density...")
    mf = dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()
    dm_scf = mf.make_rdm1()
    
    # Convert to NumPy if needed
    dm_scf = np.asarray(dm_scf, dtype=np.float64)
    
    print(f"   SCF energy: {mf.e_tot:.6f} a.u.")
    print(f"   SCF density shape: {dm_scf.shape}")
    print(f"   SCF electrons: {np.trace(dm_scf).item():.2f}")
    
    # 2. Calculate promolecule density (superposition of atomic densities)
    print("\n2. Calculating promolecule density (SAD - Superposition of Atomic Densities)...")
    
    # Use init_guess_by_atom to get atomic density superposition
    # This creates density from non-interacting atoms
    dm_atom = scf.hf.init_guess_by_atom(mol)
    dm_atom = np.asarray(dm_atom, dtype=np.float64)
    
    print(f"   Promolecule density shape: {dm_atom.shape}")
    print(f"   Promolecule electrons: {np.trace(dm_atom).item():.2f}")
    
    # 3. Calculate deformation density
    print("\n3. Calculating deformation density...")
    dm_deformation = dm_scf - dm_atom
    
    print(f"   Deformation density shape: {dm_deformation.shape}")
    print(f"   Deformation integral: {np.trace(dm_deformation).item():.6f} (should be ~0)")
    print(f"   Max deformation: {np.max(np.abs(dm_deformation)):.6f}")
    
    # 4. Save cube files
    print("\n4. Saving cube files...")
    os.makedirs('test_output', exist_ok=True)
    
    cubegen.density(mol, 'test_output/h2o_rks_cpu_scf_density.cube', dm_scf)
    print("   ✓ SCF density: test_output/h2o_rks_cpu_scf_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_rks_cpu_promolecule_density.cube', dm_atom)
    print("   ✓ Promolecule density: test_output/h2o_rks_cpu_promolecule_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_rks_cpu_deformation_density.cube', dm_deformation)
    print("   ✓ Deformation density: test_output/h2o_rks_cpu_deformation_density.cube")
    
    print("\n✓ RKS CPU test completed successfully!")
    return True

def test_deformation_density_uks_cpu():
    """Test deformation density for UKS (open-shell) on CPU"""
    print("\n" + "="*70)
    print("TEST 2: UKS (Open-shell) - CPU")
    print("="*70)
    
    # Create H2O+ cation (charge +1, open-shell)
    mol = gto.M(
        atom='O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692',
        basis='6-31g',
        charge=1,
        spin=1,  # 2S = 1, so doublet
        verbose=0
    )
    
    print(f"Molecule: H2O+ (cation)")
    print(f"Electrons: {mol.nelectron}")
    print(f"Basis functions: {mol.nao_nr()}")
    
    # 1. Calculate SCF density
    print("\n1. Calculating UKS SCF density...")
    mf = dft.UKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()
    dm_scf = mf.make_rdm1()
    
    # Handle UKS format: shape (2, nao, nao)
    if hasattr(dm_scf, 'ndim') and dm_scf.ndim == 3 and dm_scf.shape[0] == 2:
        dm_scf = np.asarray(dm_scf, dtype=np.float64)
        dm_scf_alpha = dm_scf[0]
        dm_scf_beta = dm_scf[1]
        dm_scf_total = dm_scf_alpha + dm_scf_beta
        
        print(f"   SCF energy: {mf.e_tot:.6f} a.u.")
        print(f"   SCF density shape: {dm_scf.shape}")
        print(f"   Alpha electrons: {np.trace(dm_scf_alpha).item():.2f}")
        print(f"   Beta electrons: {np.trace(dm_scf_beta).item():.2f}")
        print(f"   Total electrons: {np.trace(dm_scf_total).item():.2f}")
    else:
        raise ValueError(f"Unexpected UKS density format: {type(dm_scf)}, shape: {getattr(dm_scf, 'shape', 'N/A')}")
    
    # 2. Calculate promolecule density
    print("\n2. Calculating promolecule density (SAD)...")
    
    # For UKS, init_guess_by_atom returns shape (2, nao, nao)
    dm_atom = scf.uhf.init_guess_by_atom(mol)
    dm_atom = np.asarray(dm_atom, dtype=np.float64)
    
    if dm_atom.ndim == 3 and dm_atom.shape[0] == 2:
        dm_atom_alpha = dm_atom[0]
        dm_atom_beta = dm_atom[1]
        dm_atom_total = dm_atom_alpha + dm_atom_beta
        
        print(f"   Promolecule density shape: {dm_atom.shape}")
        print(f"   Alpha electrons: {np.trace(dm_atom_alpha).item():.2f}")
        print(f"   Beta electrons: {np.trace(dm_atom_beta).item():.2f}")
        print(f"   Total electrons: {np.trace(dm_atom_total).item():.2f}")
    else:
        raise ValueError(f"Unexpected promolecule density format")
    
    # 3. Calculate deformation density
    print("\n3. Calculating deformation density...")
    dm_deformation_total = dm_scf_total - dm_atom_total
    
    print(f"   Deformation density shape: {dm_deformation_total.shape}")
    print(f"   Deformation integral: {np.trace(dm_deformation_total).item():.6f} (should be ~0)")
    print(f"   Max deformation: {np.max(np.abs(dm_deformation_total)):.6f}")
    
    # 4. Save cube files
    print("\n4. Saving cube files...")
    
    cubegen.density(mol, 'test_output/h2o_uks_cpu_scf_density.cube', dm_scf_total)
    print("   ✓ SCF density: test_output/h2o_uks_cpu_scf_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_uks_cpu_promolecule_density.cube', dm_atom_total)
    print("   ✓ Promolecule density: test_output/h2o_uks_cpu_promolecule_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_uks_cpu_deformation_density.cube', dm_deformation_total)
    print("   ✓ Deformation density: test_output/h2o_uks_cpu_deformation_density.cube")
    
    print("\n✓ UKS CPU test completed successfully!")
    return True

def test_deformation_density_rks_gpu():
    """Test deformation density for RKS on GPU"""
    if not GPU_AVAILABLE:
        print("\n✗ Skipping GPU RKS test (GPU4PySCF not available)")
        return False
    
    print("\n" + "="*70)
    print("TEST 3: RKS (Closed-shell) - GPU")
    print("="*70)
    
    # Create H2O molecule (neutral, closed-shell)
    mol = gto.M(
        atom='O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692',
        basis='6-31g',
        charge=0,
        spin=0,
        verbose=0
    )
    
    print(f"Molecule: H2O (neutral)")
    print(f"Electrons: {mol.nelectron}")
    print(f"Basis functions: {mol.nao_nr()}")
    
    # 1. Calculate SCF density on GPU
    print("\n1. Calculating GPU SCF density...")
    mf = gpu_dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()
    dm_scf = mf.make_rdm1()
    
    # Convert CuPy to NumPy
    if hasattr(dm_scf, 'get'):
        dm_scf = dm_scf.get()
    dm_scf = np.asarray(dm_scf, dtype=np.float64)
    
    print(f"   SCF energy: {mf.e_tot:.6f} a.u.")
    print(f"   SCF density shape: {dm_scf.shape}")
    print(f"   SCF electrons: {np.trace(dm_scf).item():.2f}")
    
    # 2. Calculate promolecule density (use CPU version)
    print("\n2. Calculating promolecule density (SAD)...")
    dm_atom = scf.hf.init_guess_by_atom(mol)
    dm_atom = np.asarray(dm_atom, dtype=np.float64)
    
    print(f"   Promolecule density shape: {dm_atom.shape}")
    print(f"   Promolecule electrons: {np.trace(dm_atom).item():.2f}")
    
    # 3. Calculate deformation density
    print("\n3. Calculating deformation density...")
    dm_deformation = dm_scf - dm_atom
    
    print(f"   Deformation density shape: {dm_deformation.shape}")
    print(f"   Deformation integral: {np.trace(dm_deformation).item():.6f} (should be ~0)")
    print(f"   Max deformation: {np.max(np.abs(dm_deformation)):.6f}")
    
    # 4. Save cube files
    print("\n4. Saving cube files...")
    
    cubegen.density(mol, 'test_output/h2o_rks_gpu_scf_density.cube', dm_scf)
    print("   ✓ SCF density: test_output/h2o_rks_gpu_scf_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_rks_gpu_promolecule_density.cube', dm_atom)
    print("   ✓ Promolecule density: test_output/h2o_rks_gpu_promolecule_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_rks_gpu_deformation_density.cube', dm_deformation)
    print("   ✓ Deformation density: test_output/h2o_rks_gpu_deformation_density.cube")
    
    print("\n✓ RKS GPU test completed successfully!")
    return True

def test_deformation_density_uks_gpu():
    """Test deformation density for UKS on GPU"""
    if not GPU_AVAILABLE:
        print("\n✗ Skipping GPU UKS test (GPU4PySCF not available)")
        return False
    
    print("\n" + "="*70)
    print("TEST 4: UKS (Open-shell) - GPU")
    print("="*70)
    
    # Create H2O+ cation (charge +1, open-shell)
    mol = gto.M(
        atom='O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692',
        basis='6-31g',
        charge=1,
        spin=1,
        verbose=0
    )
    
    print(f"Molecule: H2O+ (cation)")
    print(f"Electrons: {mol.nelectron}")
    print(f"Basis functions: {mol.nao_nr()}")
    
    # 1. Calculate SCF density on GPU
    print("\n1. Calculating GPU UKS SCF density...")
    mf = gpu_dft.UKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()
    dm_scf = mf.make_rdm1()
    
    # Convert CuPy to NumPy
    if hasattr(dm_scf, 'get'):
        dm_scf = dm_scf.get()
    dm_scf = np.asarray(dm_scf, dtype=np.float64)
    
    # Handle UKS format
    if dm_scf.ndim == 3 and dm_scf.shape[0] == 2:
        dm_scf_alpha = dm_scf[0]
        dm_scf_beta = dm_scf[1]
        dm_scf_total = dm_scf_alpha + dm_scf_beta
        
        print(f"   SCF energy: {mf.e_tot:.6f} a.u.")
        print(f"   SCF density shape: {dm_scf.shape}")
        print(f"   Alpha electrons: {np.trace(dm_scf_alpha).item():.2f}")
        print(f"   Beta electrons: {np.trace(dm_scf_beta).item():.2f}")
        print(f"   Total electrons: {np.trace(dm_scf_total).item():.2f}")
    else:
        raise ValueError(f"Unexpected UKS density format")
    
    # 2. Calculate promolecule density
    print("\n2. Calculating promolecule density (SAD)...")
    dm_atom = scf.uhf.init_guess_by_atom(mol)
    dm_atom = np.asarray(dm_atom, dtype=np.float64)
    
    if dm_atom.ndim == 3 and dm_atom.shape[0] == 2:
        dm_atom_alpha = dm_atom[0]
        dm_atom_beta = dm_atom[1]
        dm_atom_total = dm_atom_alpha + dm_atom_beta
        
        print(f"   Promolecule density shape: {dm_atom.shape}")
        print(f"   Alpha electrons: {np.trace(dm_atom_alpha).item():.2f}")
        print(f"   Beta electrons: {np.trace(dm_atom_beta).item():.2f}")
        print(f"   Total electrons: {np.trace(dm_atom_total).item():.2f}")
    else:
        raise ValueError(f"Unexpected promolecule density format")
    
    # 3. Calculate deformation density
    print("\n3. Calculating deformation density...")
    dm_deformation_total = dm_scf_total - dm_atom_total
    
    print(f"   Deformation density shape: {dm_deformation_total.shape}")
    print(f"   Deformation integral: {np.trace(dm_deformation_total).item():.6f} (should be ~0)")
    print(f"   Max deformation: {np.max(np.abs(dm_deformation_total)):.6f}")
    
    # 4. Save cube files
    print("\n4. Saving cube files...")
    
    cubegen.density(mol, 'test_output/h2o_uks_gpu_scf_density.cube', dm_scf_total)
    print("   ✓ SCF density: test_output/h2o_uks_gpu_scf_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_uks_gpu_promolecule_density.cube', dm_atom_total)
    print("   ✓ Promolecule density: test_output/h2o_uks_gpu_promolecule_density.cube")
    
    cubegen.density(mol, 'test_output/h2o_uks_gpu_deformation_density.cube', dm_deformation_total)
    print("   ✓ Deformation density: test_output/h2o_uks_gpu_deformation_density.cube")
    
    print("\n✓ UKS GPU test completed successfully!")
    return True

if __name__ == '__main__':
    print("="*70)
    print("DEFORMATION DENSITY TEST SUITE")
    print("="*70)
    print("\nThis script tests calculation of deformation density:")
    print("  Deformation density = SCF density - Promolecule density")
    print("\nWhere:")
    print("  - SCF density: from converged DFT calculation (with e-e interaction)")
    print("  - Promolecule density: superposition of atomic densities (no e-e interaction)")
    print("  - Deformation density: shows charge redistribution due to bonding")
    
    results = {}
    
    # Run tests
    try:
        results['RKS_CPU'] = test_deformation_density_rks_cpu()
    except Exception as e:
        print(f"\n✗ RKS CPU test failed: {e}")
        results['RKS_CPU'] = False
    
    try:
        results['UKS_CPU'] = test_deformation_density_uks_cpu()
    except Exception as e:
        print(f"\n✗ UKS CPU test failed: {e}")
        results['UKS_CPU'] = False
    
    try:
        results['RKS_GPU'] = test_deformation_density_rks_gpu()
    except Exception as e:
        print(f"\n✗ RKS GPU test failed: {e}")
        results['RKS_GPU'] = False
    
    try:
        results['UKS_GPU'] = test_deformation_density_uks_gpu()
    except Exception as e:
        print(f"\n✗ UKS GPU test failed: {e}")
        results['UKS_GPU'] = False
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    for test_name, result in results.items():
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"{test_name:15s}: {status}")
    
    print("\nOutput files saved to: test_output/")
    print("\nVisualization:")
    print("  vmd test_output/h2o_*_deformation_density.cube")
    print("  - Positive (red): electron accumulation (bonding regions)")
    print("  - Negative (blue): electron depletion (atomic cores)")
    print("="*70)
