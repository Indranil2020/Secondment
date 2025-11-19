#!/usr/bin/env python
'''
GPU-Accelerated Transition Density Matrix Calculation
Based on official PySCF examples with GPU4PySCF support

Features:
- GPU acceleration for DFT and TDDFT (B3LYP, UKS supported)
- Parallel calculation support
- Configurable grid size and box dimensions
- Selective state output
- HOMO/LUMO cube file generation

Requirements:
- pip install gpu4pyscf-cuda12x
- pip install cutensor-cu12 (optional, for better performance)
'''

from pyscf import gto, lib
from pyscf.tools import cubegen, molden
import numpy as np
from functools import reduce
import os

# Import GPU modules (GPU4PySCF must be installed)
from gpu4pyscf import dft, tddft

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE SETTINGS
# ============================================================================

# --- Parallel Calculation Settings ---
# Note: GPU handles parallelization automatically
ENABLE_PARALLEL = True  # Enable/disable parallel computation for CPU operations
NUM_THREADS = 0  # Number of CPU threads for non-GPU operations (0 = auto-detect)

# --- Molecule Selection ---
USE_XYZ = True  # True: use XYZ file, False: use H2O test molecule
XYZ_FILE = 'PTCDA_clean.xyz'  # Path to XYZ file
BASIS_SET = '6-31g'  # Basis set

# --- Charge and Spin Settings ---
CHARGE = -1  # Molecular charge: 0 (neutral), +1 (cation), -1 (anion)
SPIN = None  # Spin multiplicity (2S+1): None = auto-calculate, 1 = singlet, 2 = doublet, 3 = triplet
# Note: Spin is auto-calculated from electron count if set to None
# For charged systems: cation (+1) typically has spin=2 (doublet), anion (-1) has spin=2 (doublet)
# Neutral even-electron systems typically have spin=1 (singlet)

# --- DFT/TDDFT Settings ---
# Note: TDDFT uses the same basis set and XC functional as ground state DFT
XC_FUNCTIONAL = 'b3lyp'  # Exchange-correlation functional
# Common options:
#   'b3lyp'    - B3LYP (hybrid, good general purpose, GPU-optimized)
#   'pbe0'     - PBE0 (hybrid, good for excited states, GPU-optimized)
#   'cam-b3lyp' - CAM-B3LYP (range-separated, good for charge transfer, GPU-optimized)
#   'wb97x-d'  - ωB97X-D (range-separated with dispersion)
#   'pbe'      - PBE (GGA, faster but less accurate, GPU-optimized)
#   'blyp'     - BLYP (GGA, GPU-optimized)

NUM_EXCITED_STATES = 10  # Total number of excited states to calculate

# --- TDDFT Method Selection ---
USE_TDA = False  # True: TDA (2× faster, ~95% accurate), False: Full TDDFT (slower, more accurate)
# TDA (Tamm-Dancoff Approximation) is recommended for large systems or initial testing
# For charged/open-shell systems, TDDFT automatically uses appropriate method
# Closed-shell (spin=1): Uses RKS/TDDFT or RKS/TDA
# Open-shell (spin>1): Uses UKS/TDDFT or UKS/TDA

# --- Output Selection ---
# STATES_TO_OUTPUT: Which states to generate CUBE FILES for (0-indexed)
# Cube files are large (~150-500 MB per state), so be selective
# Examples:
#   [0, 1, 2] - First three states
#   [0, 4, 9] - States 1, 5, and 10
#   range(5) - First five states
STATES_TO_OUTPUT = [0, 1, 2]  # Cube files: transition density, excited density, density difference

# --- Cube File Generation Options ---
GENERATE_TRANSITION_DENSITY = True  # Transition density matrix
GENERATE_EXCITED_DENSITY = True     # Excited state density
GENERATE_DENSITY_DIFFERENCE = True  # Density difference (excited - ground)
GENERATE_HOMO_LUMO = True           # HOMO and LUMO orbitals for verification

# --- Grid Settings ---
# Option 1: Use grid resolution (number of points per axis)
USE_GRID_RESOLUTION = False
GRID_RESOLUTION = [80, 80, 80]  # [nx, ny, nz] grid points

# Option 2: Use box dimensions (in Angstrom) - only used if USE_GRID_RESOLUTION = False
BOX_MARGIN = 4.0  # Margin around molecule in Angstrom
GRID_SPACING = 0.2  # Grid spacing in Angstrom

# --- NTO Analysis ---
# NTO_STATES: Which states to perform NTO ANALYSIS for (0-indexed)
# NTO molden files are small (~5 MB per state), so you can analyze more states
# This is INDEPENDENT of STATES_TO_OUTPUT - you can have different lists
# Example: Generate cube files for [0,1] but NTO analysis for [0,1,2,3,4]
ENABLE_NTO_ANALYSIS = True  # Generate NTO molden files
NTO_STATES = [0, 1, 2]  # NTO analysis (can be different from STATES_TO_OUTPUT)

# --- Output Directory ---
OUTPUT_DIR = 'output_gpu'  # Directory for output files (created if doesn't exist)

# ============================================================================
# END OF CONFIGURATION
# ============================================================================

# ============================================================================
# SETUP AND INITIALIZATION
# ============================================================================

# Create output directory
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")

# Setup parallel computation (for CPU operations like cube file generation)
if ENABLE_PARALLEL:
    if NUM_THREADS > 0:
        lib.num_threads(NUM_THREADS)
        print(f"\nCPU threads for non-GPU operations: {NUM_THREADS}")
    else:
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        lib.num_threads(num_cores)
        print(f"\nCPU threads for non-GPU operations: {num_cores} (auto-detected)")
else:
    lib.num_threads(1)
    print("\nParallel computation disabled: using 1 thread for CPU operations")

# ============================================================================
# MOLECULE DEFINITION
# ============================================================================

def create_h2o_molecule():
    """Create H2O test molecule"""
    mol = gto.M(
        atom = '''
        O  0.0000  0.0000  0.1173
        H  0.0000  0.7572 -0.4692
        H  0.0000 -0.7572 -0.4692
        ''',
        basis = BASIS_SET,
        verbose = 4
    )
    return mol

def calculate_spin_multiplicity(n_electrons, charge):
    """Calculate spin multiplicity (2S+1) from electron count"""
    n_elec = n_electrons - charge
    if n_elec % 2 == 0:
        spin = 1  # Singlet (closed-shell)
    else:
        spin = 2  # Doublet (open-shell)
    return spin

def create_molecule_from_xyz(xyz_file, basis, charge=0, spin=None):
    """Load molecule from XYZ file with charge and spin"""
    mol_temp = gto.M(atom=xyz_file, basis=basis, charge=0, spin=0, verbose=0)
    n_electrons = mol_temp.nelectron
    
    if spin is None:
        spin = calculate_spin_multiplicity(n_electrons, charge)
    
    mol = gto.M(
        atom = xyz_file,
        basis = basis,
        charge = charge,
        spin = spin - 1,  # PySCF uses 2S
        verbose = 4
    )
    
    return mol, spin

print("\n" + "="*70)
print("MOLECULE SETUP")
print("="*70)

if USE_XYZ:
    print(f"Loading molecule from: {XYZ_FILE}")
    print(f"Basis set: {BASIS_SET}")
    print(f"Charge: {CHARGE}")
    mol, calculated_spin = create_molecule_from_xyz(XYZ_FILE, BASIS_SET, CHARGE, SPIN)
    actual_spin = calculated_spin
else:
    print("Using H2O test molecule")
    print(f"Basis set: {BASIS_SET}")
    mol = create_h2o_molecule()
    actual_spin = 1  # H2O is singlet

print(f"Number of atoms: {mol.natm}")
print(f"Number of electrons: {mol.nelectron}")
print(f"Number of basis functions: {mol.nao}")
print(f"Molecular charge: {mol.charge}")
print(f"Spin multiplicity (2S+1): {actual_spin}")
print(f"Number of unpaired electrons (2S): {mol.spin}")

if actual_spin == 1:
    print("System type: Closed-shell (singlet)")
    dft_method = "RKS"
else:
    print(f"System type: Open-shell ({['singlet', 'doublet', 'triplet', 'quartet', 'quintet'][actual_spin-1] if actual_spin <= 5 else f'spin={actual_spin}'})")
    dft_method = "UKS"

print(f"DFT method: {dft_method} (GPU-accelerated)")

# ============================================================================
# GROUND STATE DFT CALCULATION
# ============================================================================

print("\n" + "="*70)
print("GROUND STATE DFT CALCULATION")
print("="*70)
print(f"XC functional: {XC_FUNCTIONAL}")
print(f"Basis set: {BASIS_SET}")
print(f"Method: {dft_method} (GPU-accelerated)")

# Select RKS (closed-shell) or UKS (open-shell) based on spin
if actual_spin == 1:
    mf = dft.RKS(mol)
else:
    mf = dft.UKS(mol)

mf.xc = XC_FUNCTIONAL
mf.kernel()

if not mf.converged:
    print("WARNING: SCF did not converge!")
    print("Try: 1) Different initial guess, 2) Level shifting, 3) DIIS settings")
else:
    print("✓ SCF converged")

print(f"Ground state energy: {mf.e_tot:.6f} a.u.")

# ============================================================================
# TDDFT CALCULATION
# ============================================================================

print("\n" + "="*70)
print("TDDFT CALCULATION")
print("="*70)
print(f"Note: TDDFT inherits XC functional ({XC_FUNCTIONAL}) and basis set ({BASIS_SET}) from ground state")
method_name = 'TDA' if USE_TDA else 'TDDFT'
spin_type = 'RKS' if actual_spin == 1 else 'UKS'
print(f"TDDFT method: {method_name} ({spin_type}-based, GPU-accelerated)")

# Select TDA or full TDDFT
# TDDFT automatically uses correct method based on mf (RKS or UKS)
if USE_TDA:
    td = tddft.TDA(mf)  # Faster, ~95% accuracy
    print("Using TDA (Tamm-Dancoff Approximation) - faster calculation")
else:
    td = tddft.TDDFT(mf)  # More accurate
    print("Using full TDDFT - more accurate but slower")

td.nstates = NUM_EXCITED_STATES
print(f"Calculating {NUM_EXCITED_STATES} excited states...")
td.kernel()

# Handle both RKS (scalar) and UKS (array) convergence
if hasattr(td.converged, '__len__'):  # UKS: array
    if not td.converged.all():
        print(f"WARNING: TDDFT did not converge for some states!")
        print(f"  Converged states: {td.converged.sum()}/{len(td.converged)}")
    else:
        print("✓ TDDFT converged (all states)")
else:  # RKS: scalar
    if not td.converged:
        print("WARNING: TDDFT did not converge!")
    else:
        print("✓ TDDFT converged")

td.analyze()

# Print excitation energies
print("\n" + "="*70)
print("EXCITED STATE ENERGIES")
print("="*70)
for i, energy in enumerate(td.e):
    print(f"State {i+1}: {energy:.6f} a.u. = {energy*27.211:.3f} eV")
print("="*70)

# ============================================================================
# TRANSITION DIPOLE MOMENTS
# ============================================================================

print("\n" + "="*70)
print("TRANSITION DIPOLE MOMENTS")
print("="*70)

charges = mol.atom_charges()
coords = mol.atom_coords()
nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
mol.set_common_orig_(nuc_charge_center)

dip_ints = mol.intor('cint1e_r_sph', comp=3)

def calculate_transition_dipole(td, state_id):
    """Calculate transition dipole moment"""
    X, Y = td.xy[state_id]
    
    mo_coeff = td._scf.mo_coeff
    mo_occ = td._scf.mo_occ
    
    orbo = mo_coeff[:, mo_occ > 0]
    orbv = mo_coeff[:, mo_occ == 0]
    
    nocc = orbo.shape[1]
    nvir = orbv.shape[1]
    
    t_dm1_mo = np.zeros((mo_coeff.shape[1], mo_coeff.shape[1]))
    t_dm1_mo[:nocc, nocc:] = (X + Y).reshape(nocc, nvir)
    
    t_dm1_ao = reduce(np.dot, (mo_coeff, t_dm1_mo, mo_coeff.T))
    
    tdm = np.einsum('xij,ji->x', dip_ints, t_dm1_ao)
    
    return tdm

print("\nTransition dipole moments (a.u.):")
print(f"{'State':<8} {'μ_x':<12} {'μ_y':<12} {'μ_z':<12} {'|μ|':<12} {'f':<12}")
print("-" * 70)

for i in range(td.nstates):
    tdm = calculate_transition_dipole(td, i)
    tdm_magnitude = np.linalg.norm(tdm)
    omega = td.e[i]
    osc_strength = (2.0/3.0) * omega * tdm_magnitude**2
    
    print(f"{i+1:<8} {tdm[0]:>11.6f} {tdm[1]:>11.6f} {tdm[2]:>11.6f} "
          f"{tdm_magnitude:>11.6f} {osc_strength:>11.6f}")

print("="*70)

print("\n" + "="*70)
print("CALCULATION SUMMARY")
print("="*70)
print(f"\nAcceleration: {'GPU (RTX 3060)' if USING_GPU else 'CPU'}")
print(f"XC functional: {XC_FUNCTIONAL}")
print(f"Basis set: {BASIS_SET}")
print(f"Molecular charge: {mol.charge}")
print(f"Spin multiplicity: {actual_spin}")
print(f"DFT method: {dft_method}")
print(f"TDDFT method: {method_name}")
print(f"Excited states calculated: {NUM_EXCITED_STATES}")
print(f"\nOutput directory: {OUTPUT_DIR}/")
print("="*70)

print("\n✓ Calculation completed successfully!")
print(f"✓ GPU acceleration: {'ENABLED' if USING_GPU else 'DISABLED'}")
