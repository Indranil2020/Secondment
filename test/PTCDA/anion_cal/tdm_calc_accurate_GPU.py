#!/usr/bin/env python
'''
GPU-Accelerated Transition Density Matrix Calculation
Based on official PySCF examples with GPU4PySCF support

Features:
- GPU acceleration for DFT and TDDFT (B3LYP, UKS fully supported)
- Parallel calculation support
- Configurable grid size and box dimensions
- Selective state output
- HOMO/LUMO cube file generation

Requirements:
- pip install gpu4pyscf-cuda12x
- pip install cutensor-cu12 (optional, for 10-20% better performance)
'''

from pyscf import gto, lib
from pyscf.tools import cubegen, molden
from gpu4pyscf import dft
from gpu4pyscf.tdscf import rks as gpu_tdrks, uks as gpu_tduks
import numpy as np
from functools import reduce
import os

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE SETTINGS
# ============================================================================

# --- Parallel Calculation Settings ---
# Note: GPU handles DFT/TDDFT parallelization automatically
ENABLE_PARALLEL = True  # Enable/disable parallel computation for CPU operations
NUM_THREADS = 0

# --- Molecule Selection ---
USE_XYZ = True
# XYZ_FILE = 'H2O.xyz'  # Path to XYZ file
XYZ_FILE = 'H2O.xyz'
BASIS_SET = '6-31g'

# --- Charge and Spin Settings ---
CHARGE = -1
SPIN = None
# Note: Spin is auto-calculated from electron count if set to None
# For charged systems: cation (+1) typically has spin=2 (doublet), anion (-1) has spin=2 (doublet)
# Neutral even-electron systems typically have spin=1 (singlet)

# --- DFT/TDDFT Settings ---
# Note: TDDFT uses the same basis set and XC functional as ground state DFT
XC_FUNCTIONAL = 'b3lyp'
# Common options:
#   'b3lyp'    - B3LYP (hybrid, good general purpose)
#   'pbe0'     - PBE0 (hybrid, good for excited states)
#   'cam-b3lyp' - CAM-B3LYP (range-separated, good for charge transfer)
#   'wb97x-d'  - ωB97X-D (range-separated with dispersion)
#   'pbe'      - PBE (GGA, faster but less accurate)
#   'blyp'     - BLYP (GGA)

NUM_EXCITED_STATES = 10

# --- TDDFT Method Selection ---
USE_TDA = False
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
STATES_TO_OUTPUT = [0, 1, 2]

# --- Cube File Generation Options ---
GENERATE_TRANSITION_DENSITY = True
GENERATE_EXCITED_DENSITY = True
GENERATE_DENSITY_DIFFERENCE = True
GENERATE_HOMO_LUMO = True

# --- Grid Settings ---
# Option 1: Use grid resolution (number of points per axis)
USE_GRID_RESOLUTION = False
GRID_RESOLUTION = [80, 80, 80]

# Option 2: Use box dimensions (in Angstrom) - only used if USE_GRID_RESOLUTION = False
BOX_MARGIN = 4.0
GRID_SPACING = 0.2

# --- NTO Analysis ---
# NTO_STATES: Which states to perform NTO ANALYSIS for (0-indexed)
# NTO molden files are small (~5 MB per state), so you can analyze more states
# This is INDEPENDENT of STATES_TO_OUTPUT - you can have different lists
# Example: Generate cube files for [0,1] but NTO analysis for [0,1,2,3,4]
ENABLE_NTO_ANALYSIS = True
NTO_STATES = [0, 1, 2]

# --- Transition Contribution Analysis ---
# Analyze which orbital pairs (i→a) contribute to each excited state
# Shows percentage contribution and generates cube files for dominant pairs
ENABLE_CONTRIBUTION_ANALYSIS = True
CONTRIBUTION_STATES = [0, 1, 2]
CONTRIBUTION_THRESHOLD = 0.01
TOP_N_CONTRIBUTIONS = 10
GENERATE_PAIR_CUBES = True
MAX_PAIRS_PER_STATE = 3
PAIR_CONTRIBUTION_CUTOFF = 0.05

# --- Ground State Density and Potential ---
# Generate cube files for ground state charge density and electrostatic potential
GENERATE_GROUND_STATE_DENSITY = True
GENERATE_ELECTROSTATIC_POTENTIAL = True

# --- Output Directory ---
OUTPUT_DIR = 'output_gpu_charge-1'

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

# Setup parallel computation
if ENABLE_PARALLEL:
    if NUM_THREADS > 0:
        lib.num_threads(NUM_THREADS)
        print(f"\nParallel computation enabled: {NUM_THREADS} threads")
    else:
        # Auto-detect number of cores
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        lib.num_threads(num_cores)
        print(f"\nParallel computation enabled: {num_cores} threads (auto-detected)")
else:
    lib.num_threads(1)
    print("\nParallel computation disabled: using 1 thread")

# ============================================================================
# 1. MOLECULE DEFINITION
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
        verbose = 0
    )
    return mol

def calculate_spin_multiplicity(n_electrons, charge):
    """
    Calculate spin multiplicity (2S+1) from electron count.
    
    Args:
        n_electrons: Total number of electrons in neutral molecule
        charge: Molecular charge
    
    Returns:
        spin: Spin multiplicity (2S+1)
    """
    # Adjust electron count for charge
    n_elec = n_electrons - charge
    
    # For even number of electrons: singlet (spin=1)
    # For odd number of electrons: doublet (spin=2)
    if n_elec % 2 == 0:
        spin = 1  # Singlet (closed-shell)
    else:
        spin = 2  # Doublet (open-shell)
    
    return spin

def create_molecule_from_xyz(xyz_file, basis, charge=0, spin=None):
    """
    Load molecule from XYZ file with charge and spin.
    
    Args:
        xyz_file: Path to XYZ file
        basis: Basis set
        charge: Molecular charge
        spin: Spin multiplicity (None = auto-calculate)
    """
    # First, create molecule to get electron count
    mol_temp = gto.M(atom=xyz_file, basis=basis, charge=0, spin=0, verbose=0)
    n_electrons = mol_temp.nelectron
    
    # Calculate spin if not provided
    if spin is None:
        spin = calculate_spin_multiplicity(n_electrons, charge)
    
    # Create final molecule with charge and spin
    mol = gto.M(
        atom = xyz_file,
        basis = basis,
        charge = charge,
        spin = spin - 1,  # PySCF uses 2S (number of unpaired electrons), not 2S+1
        verbose = 0
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
# 2. GROUND STATE DFT CALCULATION
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
# 2A. GROUND STATE DENSITY AND ELECTROSTATIC POTENTIAL
# ============================================================================

if GENERATE_GROUND_STATE_DENSITY or GENERATE_ELECTROSTATIC_POTENTIAL:
    print("\n" + "="*70)
    print("GROUND STATE DENSITY AND POTENTIAL")
    print("="*70)
    
    # Calculate ground state density matrix
    dm = mf.make_rdm1()
    
    # Convert CuPy to NumPy first (GPU4PySCF returns CuPy arrays)
    if hasattr(dm, 'get'):
        dm = dm.get()
    
    # Handle UKS (sum alpha and beta densities)
    # GPU4PySCF UKS returns shape (2, nao, nao) or tuple of (dm_alpha, dm_beta)
    # CPU PySCF UKS returns tuple of (dm_alpha, dm_beta)
    if isinstance(dm, tuple):
        # Standard PySCF UKS: tuple of (dm_alpha, dm_beta)
        dm_alpha, dm_beta = dm
        dm_alpha = np.asarray(dm_alpha, dtype=np.float64)
        dm_beta = np.asarray(dm_beta, dtype=np.float64)
        
        dm_total = dm_alpha + dm_beta
        
        n_alpha = np.trace(dm_alpha).item()  # .item() converts to Python scalar
        n_beta = np.trace(dm_beta).item()
        
        print("System type: UKS (open-shell)")
        print(f"  Alpha electrons: {n_alpha:.2f}")
        print(f"  Beta electrons: {n_beta:.2f}")
    elif dm.ndim == 3 and dm.shape[0] == 2:
        # GPU4PySCF UKS: array with shape (2, nao, nao)
        dm = np.asarray(dm, dtype=np.float64)
        dm_alpha = dm[0]
        dm_beta = dm[1]
        
        dm_total = dm_alpha + dm_beta
        
        n_alpha = np.trace(dm_alpha).item()
        n_beta = np.trace(dm_beta).item()
        
        print("System type: UKS (open-shell)")
        print(f"  Alpha electrons: {n_alpha:.2f}")
        print(f"  Beta electrons: {n_beta:.2f}")
    else:
        # RKS: single 2D matrix
        dm_total = np.asarray(dm, dtype=np.float64)
        print("System type: RKS (closed-shell)")
    
    # Calculate total electrons
    total_electrons = np.trace(dm_total).item()  # .item() for scalar
    print(f"Total electrons: {total_electrons:.2f}")
    print(f"Molecular charge: {CHARGE}")
    print(f"Expected electrons: {mol.nelectron}")
    
    # Verify density matrix is proper NumPy array
    if not isinstance(dm_total, np.ndarray):
        raise TypeError(f"Density matrix is not NumPy array: {type(dm_total)}")
    
    nao = mol.nao_nr()
    if dm_total.shape != (nao, nao):
        raise ValueError(f"Density matrix shape {dm_total.shape} doesn't match AO basis {nao}x{nao}")
    
    # Generate ground state charge density cube file
    if GENERATE_GROUND_STATE_DENSITY:
        density_file = os.path.join(OUTPUT_DIR, 'ground_state_density.cube')
        print(f"\nGenerating ground state charge density...")
        try:
            cubegen.density(mol, density_file, dm_total)
            print(f"  ✓ Ground state density: {density_file}")
            print(f"    Use this to visualize total electron distribution")
        except Exception as e:
            print(f"  ✗ Failed to generate density cube: {str(e)}")
    
    # Generate electrostatic potential cube file
    if GENERATE_ELECTROSTATIC_POTENTIAL:
        esp_file = os.path.join(OUTPUT_DIR, 'electrostatic_potential.cube')
        print(f"\nGenerating electrostatic potential (ESP)...")
        try:
            # ESP = Nuclear potential + Electronic potential
            # cubegen.mep calculates the molecular electrostatic potential
            from pyscf.tools import cubegen
            cubegen.mep(mol, esp_file, dm_total)
            print(f"  ✓ Electrostatic potential: {esp_file}")
            print(f"    Use this to identify:")
            print(f"      - Nucleophilic sites (negative ESP, red)")
            print(f"      - Electrophilic sites (positive ESP, blue)")
            print(f"      - Reaction sites and molecular recognition")
        except Exception as e:
            print(f"  ✗ Failed to generate ESP cube: {str(e)}")
    
    print("="*70)

# ============================================================================
# 3. TDDFT CALCULATION
# ============================================================================

print("\n" + "="*70)
print("TDDFT CALCULATION")
print("="*70)
print(f"Note: TDDFT inherits XC functional ({XC_FUNCTIONAL}) and basis set ({BASIS_SET}) from ground state")
method_name = 'TDA' if USE_TDA else 'TDDFT'
spin_type = 'RKS' if actual_spin == 1 else 'UKS'
print(f"TDDFT method: {method_name} ({spin_type}-based, GPU-accelerated)")

# Select TDA or full TDDFT
# Use gpu4pyscf.tdscf module (RKS or UKS based on spin)
if actual_spin == 1:
    # Closed-shell: use gpu4pyscf.tdscf.rks
    if USE_TDA:
        td = gpu_tdrks.TDA(mf)  # Faster, ~95% accuracy
        print("Using TDA (Tamm-Dancoff Approximation) - faster calculation")
    else:
        td = gpu_tdrks.TDDFT(mf)  # More accurate
        print("Using full TDDFT - more accurate but slower")
else:
    # Open-shell: use gpu4pyscf.tdscf.uks
    if USE_TDA:
        td = gpu_tduks.TDA(mf)  # Faster, ~95% accuracy
        print("Using TDA (Tamm-Dancoff Approximation) - faster calculation")
    else:
        td = gpu_tduks.TDDFT(mf)  # More accurate
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

td.analyze()  # Print detailed analysis

# Print excitation energies
print("\n" + "="*70)
print("EXCITED STATE ENERGIES")
print("="*70)
for i, energy in enumerate(td.e):
    print(f"State {i+1}: {energy:.6f} a.u. = {energy*27.211:.3f} eV")
print("="*70)

# ============================================================================
# 4. TRANSITION DIPOLE MOMENTS
# ============================================================================

print("\n" + "="*70)
print("TRANSITION DIPOLE MOMENTS")
print("="*70)

# Set gauge origin to nuclear charge center
charges = mol.atom_charges()
coords = mol.atom_coords()  # in a.u.
nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
mol.set_common_orig_(nuc_charge_center)

# Calculate dipole integrals
dip_ints = mol.intor('cint1e_r_sph', comp=3)  # x, y, z components

def calculate_transition_dipole(td, state_id):
    """
    Calculate transition dipole moment for a given excited state.
    Based on PySCF example: examples/1-advanced/030-transition_dipole.py
    
    Handles both RKS (closed-shell) and UKS (open-shell) cases.
    For GPU4PySCF, converts CuPy arrays to NumPy arrays.
    
    Parameters:
    -----------
    td : TDDFT object
    state_id : int (0-indexed)
    
    Returns:
    --------
    tdm : ndarray, shape (3,)
        Transition dipole moment [x, y, z] in a.u.
    """
    # Get TDDFT amplitudes (X and Y vectors)
    X, Y = td.xy[state_id]
    
    # Get MO coefficients and occupations
    mo_coeff = td._scf.mo_coeff
    mo_occ = td._scf.mo_occ
    
    # Check if UKS by checking if X is a tuple (more reliable than mo_coeff)
    is_uks = isinstance(X, tuple)
    
    if is_uks:
        # UKS: mo_coeff and mo_occ are tuples (alpha, beta)
        # For UKS, X and Y are also tuples: ((Xa, Xb), (Ya, Yb))
        mo_coeff_a, mo_coeff_b = mo_coeff
        mo_occ_a, mo_occ_b = mo_occ
        Xa, Xb = X
        Ya, Yb = Y
        
        # Convert CuPy to NumPy if needed
        if hasattr(mo_coeff_a, 'get'):
            mo_coeff_a = mo_coeff_a.get()
            mo_coeff_b = mo_coeff_b.get()
            mo_occ_a = mo_occ_a.get()
            mo_occ_b = mo_occ_b.get()
        if hasattr(Xa, 'get'):
            Xa = Xa.get()
            Xb = Xb.get()
            Ya = Ya.get()
            Yb = Yb.get()
        
        # Xa and Xb are already separated, just need to get dimensions
        nocc_a = Xa.shape[0]
        nvir_a = Xa.shape[1]
        nocc_b = Xb.shape[0]
        nvir_b = Xb.shape[1]
        nmo_a = mo_coeff_a.shape[1]
        nmo_b = mo_coeff_b.shape[1]
        
        # Transition density matrices for alpha and beta
        t_dm1_mo_a = np.zeros((nmo_a, nmo_a))
        t_dm1_mo_a[:nocc_a, nocc_a:] = Xa + Ya
        t_dm1_ao_a = reduce(np.dot, (mo_coeff_a, t_dm1_mo_a, mo_coeff_a.T))
        
        t_dm1_mo_b = np.zeros((nmo_b, nmo_b))
        t_dm1_mo_b[:nocc_b, nocc_b:] = Xb + Yb
        t_dm1_ao_b = reduce(np.dot, (mo_coeff_b, t_dm1_mo_b, mo_coeff_b.T))
        
        # Total transition density (alpha + beta)
        t_dm1_ao = t_dm1_ao_a + t_dm1_ao_b
        
    else:
        # RKS: mo_coeff and mo_occ are arrays
        # Convert CuPy to NumPy if needed
        if hasattr(mo_coeff, 'get'):
            mo_coeff = mo_coeff.get()
            mo_occ = mo_occ.get()
        if hasattr(X, 'get'):
            X = X.get()
            Y = Y.get()
        
        orbo = mo_coeff[:, mo_occ > 0]
        orbv = mo_coeff[:, mo_occ == 0]
        nocc = orbo.shape[1]
        nvir = orbv.shape[1]
        
        # Transition density in MO basis
        t_dm1_mo = np.zeros((mo_coeff.shape[1], mo_coeff.shape[1]))
        t_dm1_mo[:nocc, nocc:] = (X + Y).reshape(nocc, nvir)
        
        # Transform to AO basis
        t_dm1_ao = reduce(np.dot, (mo_coeff, t_dm1_mo, mo_coeff.T))
    
    # Calculate transition dipole: μ = Tr(μ_op * T)
    tdm = np.einsum('xij,ji->x', dip_ints, t_dm1_ao)
    
    return tdm

# Calculate and print transition dipoles for all states
print("\nTransition dipole moments (a.u.):")
print(f"{'State':<8} {'μ_x':<12} {'μ_y':<12} {'μ_z':<12} {'|μ|':<12} {'f':<12}")
print("-" * 70)

for i in range(td.nstates):
    tdm = calculate_transition_dipole(td, i)
    tdm_magnitude = np.linalg.norm(tdm)
    
    # Oscillator strength: f = (2/3) * ω * |μ|^2
    # where ω is excitation energy in a.u.
    omega = td.e[i]
    osc_strength = (2.0/3.0) * omega * tdm_magnitude**2
    
    print(f"{i+1:<8} {tdm[0]:>11.6f} {tdm[1]:>11.6f} {tdm[2]:>11.6f} "
          f"{tdm_magnitude:>11.6f} {osc_strength:>11.6f}")

print("="*70)

# ============================================================================
# 5. TRANSITION DENSITY MATRICES
# ============================================================================

def calculate_transition_density_matrix(td, state_id):
    """
    Calculate transition density matrix between ground and excited state.
    This is the proper transition density for visualization.
    Handles both RKS and UKS, converts CuPy to NumPy.
    
    Parameters:
    -----------
    td : TDDFT object
    state_id : int (0-indexed)
    
    Returns:
    --------
    t_dm1_ao : ndarray
        Transition density matrix in AO basis
    """
    X, Y = td.xy[state_id]
    mo_coeff = td._scf.mo_coeff
    mo_occ = td._scf.mo_occ
    is_uks = isinstance(X, tuple)
    
    if is_uks:
        mo_coeff_a, mo_coeff_b = mo_coeff
        mo_occ_a, mo_occ_b = mo_occ
        Xa, Xb = X
        Ya, Yb = Y
        
        # Convert CuPy to NumPy
        if hasattr(mo_coeff_a, 'get'):
            mo_coeff_a = mo_coeff_a.get()
            mo_coeff_b = mo_coeff_b.get()
            mo_occ_a = mo_occ_a.get()
            mo_occ_b = mo_occ_b.get()
        if hasattr(Xa, 'get'):
            Xa = Xa.get()
            Xb = Xb.get()
            Ya = Ya.get()
            Yb = Yb.get()
        
        nocc_a = Xa.shape[0]
        nvir_a = Xa.shape[1]
        nocc_b = Xb.shape[0]
        nvir_b = Xb.shape[1]
        nmo_a = mo_coeff_a.shape[1]
        nmo_b = mo_coeff_b.shape[1]
        
        # Alpha
        t_dm1_mo_a = np.zeros((nmo_a, nmo_a))
        t_dm1_mo_a[:nocc_a, nocc_a:] = Xa + Ya
        t_dm1_mo_a[nocc_a:, :nocc_a] = (Xa + Ya).T
        t_dm1_ao_a = reduce(np.dot, (mo_coeff_a, t_dm1_mo_a, mo_coeff_a.T))
        
        # Beta
        t_dm1_mo_b = np.zeros((nmo_b, nmo_b))
        t_dm1_mo_b[:nocc_b, nocc_b:] = Xb + Yb
        t_dm1_mo_b[nocc_b:, :nocc_b] = (Xb + Yb).T
        t_dm1_ao_b = reduce(np.dot, (mo_coeff_b, t_dm1_mo_b, mo_coeff_b.T))
        
        t_dm1_ao = t_dm1_ao_a + t_dm1_ao_b
    else:
        # Convert CuPy to NumPy
        if hasattr(mo_coeff, 'get'):
            mo_coeff = mo_coeff.get()
            mo_occ = mo_occ.get()
        if hasattr(X, 'get'):
            X = X.get()
            Y = Y.get()
        
        orbo = mo_coeff[:, mo_occ > 0]
        orbv = mo_coeff[:, mo_occ == 0]
        nocc = orbo.shape[1]
        nvir = orbv.shape[1]
        
        t_dm1_mo = np.zeros((mo_coeff.shape[1], mo_coeff.shape[1]))
        t_dm1_mo[:nocc, nocc:] = (X + Y).reshape(nocc, nvir)
        t_dm1_mo[nocc:, :nocc] = (X + Y).reshape(nocc, nvir).T
        t_dm1_ao = reduce(np.dot, (mo_coeff, t_dm1_mo, mo_coeff.T))
    
    return t_dm1_ao

def calculate_excited_state_density(td, state_id):
    """
    Calculate excited state density matrix.
    Based on PySCF examples/tddft/22-density.py
    Handles both RKS and UKS, converts CuPy to NumPy.
    
    Parameters:
    -----------
    td : TDDFT object
    state_id : int (0-indexed)
    
    Returns:
    --------
    dm_excited : ndarray
        Excited state density matrix in AO basis
    """
    X, Y = td.xy[state_id]
    mf = td._scf
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    is_uks = isinstance(X, tuple)
    
    if is_uks:
        # UKS case - simplified approach: use transition density
        # For visualization, transition density is more meaningful
        return calculate_transition_density_matrix(td, state_id)
    else:
        # RKS case
        # Convert CuPy to NumPy
        if hasattr(mo_coeff, 'get'):
            mo_coeff = mo_coeff.get()
            mo_occ = mo_occ.get()
        if hasattr(X, 'get'):
            X = X.get()
            Y = Y.get()
        
        nocc = X.shape[0]
        
        # Density matrix changes in MO basis
        dm_oo = -np.einsum('ia,ka->ik', X.conj(), X)
        dm_oo -= np.einsum('ia,ka->ik', Y.conj(), Y)
        
        dm_vv = np.einsum('ia,ic->ac', X, X.conj())
        dm_vv += np.einsum('ia,ic->ac', Y, Y.conj())
        
        # Start with ground state density in MO basis
        dm = np.diag(mo_occ)
        
        # Add TDDFT contribution
        dm[:nocc, :nocc] += dm_oo * 2
        dm[nocc:, nocc:] += dm_vv * 2
        
        # Transform to AO basis
        dm_excited = np.einsum('pi,ij,qj->pq', mo_coeff, dm, mo_coeff.conj())
    
    return dm_excited

# ============================================================================
# 5A. TRANSITION CONTRIBUTION ANALYSIS FUNCTIONS
# ============================================================================

def get_orbital_labels(mf):
    """Get orbital labels (HOMO-n, LUMO+n) for RKS and UKS."""
    mo_occ = mf.mo_occ
    mo_coeff = mf.mo_coeff
    
    # Handle UKS (use alpha spin) and CuPy arrays
    if isinstance(mo_occ, tuple):
        mo_occ = mo_occ[0]
    if isinstance(mo_coeff, tuple):
        mo_coeff = mo_coeff[0]
    if hasattr(mo_occ, 'get'):
        mo_occ = mo_occ.get()
    if hasattr(mo_coeff, 'get'):
        mo_coeff = mo_coeff.get()
    
    # Get number of orbitals from mo_coeff shape
    nmo = mo_coeff.shape[1]
    homo_idx = np.where(mo_occ > 0)[0][-1]
    
    labels = []
    for i in range(nmo):
        if i <= homo_idx:
            offset = homo_idx - i
            labels.append('HOMO' if offset == 0 else f'HOMO-{offset}')
        else:
            offset = i - homo_idx - 1
            labels.append('LUMO' if offset == 0 else f'LUMO+{offset}')
    
    return labels, homo_idx

def analyze_transition_contributions(td, state_id, mf, threshold=0.01, top_n=10):
    """
    Analyze orbital pair contributions to a specific excited state.
    Handles both RKS and UKS systems, and CuPy arrays.
    
    Returns:
        contributions: list of (occ_idx, vir_idx, weight, label) sorted by weight
        total_weight: sum of all weights
    """
    X, Y = td.xy[state_id]
    
    # Handle UKS: X and Y are tuples (Xa, Xb), (Ya, Yb)
    # For simplicity, analyze alpha spin (dominant for most cases)
    if isinstance(X, tuple):
        X, _ = X
        Y, _ = Y
    
    # Convert CuPy to NumPy if needed
    if hasattr(X, 'get'):
        X = X.get()
    if hasattr(Y, 'get'):
        Y = Y.get()
    
    nocc, nvir = X.shape
    
    # For full TDDFT, the contribution is from (X + Y)
    amplitudes = X + Y
    
    # Calculate weights (squared amplitudes)
    weights = amplitudes ** 2
    
    # Get orbital labels
    labels, homo_idx = get_orbital_labels(mf)
    
    # Collect all contributions
    contributions = []
    for i in range(nocc):
        for a in range(nvir):
            weight = weights[i, a]
            if weight > threshold:
                occ_idx = i
                vir_idx = nocc + a
                occ_label = labels[occ_idx]
                vir_label = labels[vir_idx]
                transition_label = f"{occ_label} → {vir_label}"
                contributions.append((occ_idx, vir_idx, weight, transition_label))
    
    # Sort by weight (descending)
    contributions.sort(key=lambda x: x[2], reverse=True)
    
    # Normalize weights to sum to 1
    total_weight = sum(c[2] for c in contributions)
    if total_weight > 0:
        contributions = [(occ, vir, w/total_weight, label) 
                         for occ, vir, w, label in contributions]
    
    # Return top N
    return contributions[:top_n], total_weight

def calculate_pair_transition_density(mf, occ_idx, vir_idx):
    """
    Calculate transition density matrix for a single orbital pair i→a.
    Handles both RKS and UKS (uses alpha spin for UKS), and CuPy arrays.
    
    T_μν = C_μ^i × C_ν^a + C_μ^a × C_ν^i
    
    Returns NumPy array in AO basis.
    """
    mo_coeff = mf.mo_coeff
    
    # Handle UKS (use alpha spin)
    if isinstance(mo_coeff, tuple):
        mo_coeff_alpha = mo_coeff[0]
    else:
        mo_coeff_alpha = mo_coeff
    
    # Convert CuPy to NumPy if needed
    if hasattr(mo_coeff_alpha, 'get'):
        mo_coeff_alpha = mo_coeff_alpha.get()
    
    # Ensure we have NumPy array
    mo_coeff_alpha = np.asarray(mo_coeff_alpha)
    
    # Extract specific orbitals
    occ_mo = mo_coeff_alpha[:, occ_idx]
    vir_mo = mo_coeff_alpha[:, vir_idx]
    
    # Transition density matrix for this pair in AO basis
    # T_μν = C_μ^i × C_ν^a + C_μ^a × C_ν^i (symmetric)
    t_dm = np.outer(occ_mo, vir_mo) + np.outer(vir_mo, occ_mo)
    
    # Ensure it's a proper NumPy array with correct dtype
    t_dm = np.asarray(t_dm, dtype=np.float64)
    
    return t_dm

# ============================================================================
# 6. NATURAL TRANSITION ORBITALS (NTO) ANALYSIS
# ============================================================================

if ENABLE_NTO_ANALYSIS:
    print("\n" + "="*70)
    print("NATURAL TRANSITION ORBITAL ANALYSIS")
    print("="*70)
    
    # Filter valid NTO states
    valid_nto_states = [s for s in NTO_STATES if s < td.nstates]
    
    if not valid_nto_states:
        print("No valid NTO states selected.")
    else:
        for i in valid_nto_states:
            print(f"\nState {i+1} ({td.e[i]*27.211:.3f} eV):")
            weights, nto_coeff = td.get_nto(state=i+1, verbose=0)
            
            # Save NTO orbitals to molden format for visualization
            # For UKS, nto_coeff is a tuple (alpha, beta)
            if isinstance(nto_coeff, tuple):
                # Save alpha NTOs
                molden_file_a = os.path.join(OUTPUT_DIR, f'nto_state_{i+1}_alpha.molden')
                molden.from_mo(mol, molden_file_a, nto_coeff[0])
                print(f"  Alpha NTO orbitals saved to: {molden_file_a}")
                
                # Save beta NTOs
                molden_file_b = os.path.join(OUTPUT_DIR, f'nto_state_{i+1}_beta.molden')
                molden.from_mo(mol, molden_file_b, nto_coeff[1])
                print(f"  Beta NTO orbitals saved to: {molden_file_b}")
            else:
                # RKS case
                molden_file = os.path.join(OUTPUT_DIR, f'nto_state_{i+1}.molden')
                molden.from_mo(mol, molden_file, nto_coeff)
                print(f"  NTO orbitals saved to: {molden_file}")
    
    print("="*70)
else:
    print("\nNTO analysis disabled.")

# ============================================================================
# 6A. TRANSITION CONTRIBUTION ANALYSIS
# ============================================================================

if ENABLE_CONTRIBUTION_ANALYSIS:
    print("\n" + "="*70)
    print("TRANSITION CONTRIBUTION ANALYSIS")
    print("="*70)
    print("Analyzing orbital pair contributions to excited states...")
    
    # Filter valid states
    valid_contrib_states = [s for s in CONTRIBUTION_STATES if s < td.nstates]
    
    if not valid_contrib_states:
        print("No valid states selected for contribution analysis.")
    else:
        # Store all contributions for summary
        all_contributions = {}
        
        for state_id in valid_contrib_states:
            excitation_energy = td.e[state_id] * 27.211  # Convert to eV
            
            # Analyze contributions
            contributions, total_weight = analyze_transition_contributions(
                td, state_id, mf,
                threshold=CONTRIBUTION_THRESHOLD,
                top_n=TOP_N_CONTRIBUTIONS
            )
            
            all_contributions[state_id] = (contributions, total_weight)
            
            # Print contribution table
            print(f"\n{'='*70}")
            print(f"STATE {state_id+1}: {excitation_energy:.4f} eV")
            print(f"{'='*70}")
            print(f"{'Rank':<6} {'Transition':<20} {'Weight':<12} {'Percentage':<12} {'Cumulative':<12}")
            print(f"{'-'*70}")
            
            cumulative = 0.0
            for rank, (occ_idx, vir_idx, weight, label) in enumerate(contributions, 1):
                cumulative += weight
                print(f"{rank:<6} {label:<20} {weight:<12.6f} {weight*100:<12.2f}% {cumulative*100:<12.2f}%")
            
            print(f"{'-'*70}")
            print(f"Total weight analyzed: {total_weight:.6f}")
        
        # Save contribution tables to file
        table_file = os.path.join(OUTPUT_DIR, 'contribution_tables.txt')
        with open(table_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("ORBITAL PAIR CONTRIBUTIONS TO EXCITED STATES\n")
            f.write("="*70 + "\n\n")
            
            for state_id in valid_contrib_states:
                contributions, total_weight = all_contributions[state_id]
                excitation_energy = td.e[state_id] * 27.211
                
                f.write(f"\n{'='*70}\n")
                f.write(f"STATE {state_id+1}: {excitation_energy:.4f} eV\n")
                f.write(f"{'='*70}\n")
                f.write(f"{'Rank':<6} {'Transition':<20} {'Weight':<12} {'Percentage':<12} {'Cumulative':<12}\n")
                f.write(f"{'-'*70}\n")
                
                cumulative = 0.0
                for rank, (occ_idx, vir_idx, weight, label) in enumerate(contributions, 1):
                    cumulative += weight
                    f.write(f"{rank:<6} {label:<20} {weight:<12.6f} {weight*100:<12.2f}% {cumulative*100:<12.2f}%\n")
                
                f.write(f"{'-'*70}\n")
                f.write(f"Total weight analyzed: {total_weight:.6f}\n")
                f.write(f"{'='*70}\n\n")
        
        print(f"\n✓ Contribution tables saved to: {table_file}")
        
        # Generate cube files for dominant orbital pairs
        if GENERATE_PAIR_CUBES:
            print(f"\n{'='*70}")
            print("GENERATING ORBITAL PAIR TRANSITION DENSITY CUBE FILES")
            print(f"{'='*70}")
            print("Note: Cube files will be generated after grid parameters are calculated")
            print("      (see CUBE FILE GENERATION section below)")
    
    print("="*70)
else:
    print("\nTransition contribution analysis disabled.")

# ============================================================================
# 7. CALCULATE GRID PARAMETERS
# ============================================================================

def calculate_grid_parameters(mol, use_resolution=True, resolution=None, 
                             box_margin=4.0, grid_spacing=0.2):
    """
    Calculate grid parameters for cube file generation.
    
    Parameters:
    -----------
    mol : Mole object
    use_resolution : bool
        If True, use fixed resolution. If False, calculate from box dimensions.
    resolution : list [nx, ny, nz]
        Grid resolution (number of points per axis)
    box_margin : float
        Margin around molecule in Angstrom
    grid_spacing : float
        Grid spacing in Angstrom
    
    Returns:
    --------
    nx, ny, nz : int
        Grid resolution
    box_info : dict
        Box dimension information
    """
    # Get molecular coordinates in Angstrom
    coords = mol.atom_coords() * 0.529177  # Bohr to Angstrom
    
    # Calculate bounding box
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    mol_size = max_coords - min_coords
    
    box_info = {
        'min_coords': min_coords,
        'max_coords': max_coords,
        'mol_size': mol_size,
        'box_margin': box_margin
    }
    
    if use_resolution:
        nx, ny, nz = resolution
        actual_spacing = mol_size / np.array([nx, ny, nz])
        box_info['grid_spacing'] = actual_spacing
        box_info['total_points'] = nx * ny * nz
    else:
        # Calculate grid points from spacing and box size
        box_size = mol_size + 2 * box_margin
        nx = int(np.ceil(box_size[0] / grid_spacing))
        ny = int(np.ceil(box_size[1] / grid_spacing))
        nz = int(np.ceil(box_size[2] / grid_spacing))
        box_info['grid_spacing'] = [grid_spacing] * 3
        box_info['box_size'] = box_size
        box_info['total_points'] = nx * ny * nz
    
    return nx, ny, nz, box_info

# Calculate grid parameters
print("\n" + "="*70)
print("GRID PARAMETERS")
print("="*70)

if USE_GRID_RESOLUTION:
    nx, ny, nz, box_info = calculate_grid_parameters(
        mol, use_resolution=True, resolution=GRID_RESOLUTION
    )
    print(f"Using fixed grid resolution: {nx} × {ny} × {nz}")
    print(f"Total grid points: {box_info['total_points']:,}")
    print(f"Molecule size: {box_info['mol_size']} Å")
    print(f"Effective grid spacing: {box_info['grid_spacing']} Å")
else:
    nx, ny, nz, box_info = calculate_grid_parameters(
        mol, use_resolution=False, box_margin=BOX_MARGIN, 
        grid_spacing=GRID_SPACING
    )
    print(f"Using box dimensions with margin: {BOX_MARGIN} Å")
    print(f"Grid spacing: {GRID_SPACING} Å")
    print(f"Calculated grid resolution: {nx} × {ny} × {nz}")
    print(f"Total grid points: {box_info['total_points']:,}")
    print(f"Box size: {box_info['box_size']} Å")

print("="*70)

# ============================================================================
# 8. GENERATE HOMO/LUMO CUBE FILES
# ============================================================================

if GENERATE_HOMO_LUMO:
    print("\n" + "="*70)
    print("GENERATING HOMO/LUMO CUBE FILES")
    print("="*70)
    
    # Get HOMO and LUMO indices (handle both RKS and UKS)
    if actual_spin == 1:  # RKS (closed-shell)
        mo_occ = mf.mo_occ
        mo_coeff = mf.mo_coeff
        mo_energy = mf.mo_energy
    else:  # UKS (open-shell) - use alpha orbitals
        mo_occ = mf.mo_occ[0]  # Alpha occupation
        mo_coeff = mf.mo_coeff[0]  # Alpha coefficients
        mo_energy = mf.mo_energy[0]  # Alpha energies
        print("Note: Using alpha orbitals for HOMO/LUMO (open-shell system)")
    
    # Convert CuPy to NumPy if needed
    if hasattr(mo_occ, 'get'):
        mo_occ = mo_occ.get()
        mo_coeff = mo_coeff.get()
        mo_energy = mo_energy.get()
    
    homo_idx = np.where(mo_occ > 0)[0][-1]
    lumo_idx = np.where(mo_occ == 0)[0][0]
    
    print(f"HOMO index: {homo_idx}")
    print(f"LUMO index: {lumo_idx}")
    print(f"HOMO energy: {mo_energy[homo_idx]*27.211:.3f} eV")
    print(f"LUMO energy: {mo_energy[lumo_idx]*27.211:.3f} eV")
    print(f"HOMO-LUMO gap: {(mo_energy[lumo_idx] - mo_energy[homo_idx])*27.211:.3f} eV")
    
    # Generate HOMO cube file
    homo_file = os.path.join(OUTPUT_DIR, 'HOMO.cube')
    cubegen.orbital(mol, homo_file, mo_coeff[:, homo_idx], nx=nx, ny=ny, nz=nz)
    print(f"\n  ✓ HOMO orbital: {homo_file}")
    
    # Generate LUMO cube file
    lumo_file = os.path.join(OUTPUT_DIR, 'LUMO.cube')
    cubegen.orbital(mol, lumo_file, mo_coeff[:, lumo_idx], nx=nx, ny=ny, nz=nz)
    print(f"  ✓ LUMO orbital: {lumo_file}")
    
    # Generate HOMO-1 and LUMO+1 for additional verification
    if homo_idx > 0:
        homo1_file = os.path.join(OUTPUT_DIR, 'HOMO-1.cube')
        cubegen.orbital(mol, homo1_file, mo_coeff[:, homo_idx-1], nx=nx, ny=ny, nz=nz)
        print(f"  ✓ HOMO-1 orbital: {homo1_file}")
    
    if lumo_idx < len(mo_occ) - 1:
        lumo1_file = os.path.join(OUTPUT_DIR, 'LUMO+1.cube')
        cubegen.orbital(mol, lumo1_file, mo_coeff[:, lumo_idx+1], nx=nx, ny=ny, nz=nz)
        print(f"  ✓ LUMO+1 orbital: {lumo1_file}")
    
    print("\nVerification tip:")
    print("For the first excited state (S1), check if the transition density")
    print("resembles a HOMO→LUMO transition by comparing:")
    print("  - transition_density_state1.cube")
    print("  - HOMO.cube (electron depletion)")
    print("  - LUMO.cube (electron accumulation)")
    
    # Analytical verification: Calculate HOMO→LUMO transition density
    print("\n" + "-"*70)
    print("ANALYTICAL VERIFICATION: HOMO→LUMO Transition Density")
    print("-"*70)
    
    # Construct approximate HOMO→LUMO transition density matrix
    homo_mo = mo_coeff[:, homo_idx]
    lumo_mo = mo_coeff[:, lumo_idx]
    
    # T_approx = |HOMO⟩⟨LUMO| + |LUMO⟩⟨HOMO|
    # In AO basis: T_μν = C_μ^HOMO * C_ν^LUMO + C_μ^LUMO * C_ν^HOMO
    t_homo_lumo = np.outer(homo_mo, lumo_mo) + np.outer(lumo_mo, homo_mo)
    
    # Generate cube file for HOMO→LUMO transition density
    homo_lumo_file = os.path.join(OUTPUT_DIR, 'transition_HOMO_LUMO_analytical.cube')
    cubegen.density(mol, homo_lumo_file, t_homo_lumo, nx=nx, ny=ny, nz=nz)
    print(f"  ✓ Analytical HOMO→LUMO transition: {homo_lumo_file}")
    
    print("\nTo verify S1 is a HOMO→LUMO transition, compare:")
    print(f"  1. {os.path.join(OUTPUT_DIR, 'transition_density_state1.cube')}")
    print(f"  2. {homo_lumo_file}")
    print("\nThey should be very similar if S1 is dominated by HOMO→LUMO excitation.")
    print("You can calculate the overlap/similarity in VMD or by visual inspection.")
    
    print("="*70)
else:
    print("\nHOMO/LUMO generation disabled.")

# ============================================================================
# 9. GENERATE CUBE FILES FOR SELECTED EXCITED STATES
# ============================================================================

print("\n" + "="*70)
print("GENERATING EXCITED STATE CUBE FILES")
print("="*70)

# Filter valid states
valid_states = [s for s in STATES_TO_OUTPUT if s < td.nstates]

if not valid_states:
    print("No valid states selected for cube file generation.")
else:
    print(f"Generating cube files for states: {[s+1 for s in valid_states]}")
    print(f"Grid resolution: {nx} × {ny} × {nz}")
    
    for state_id in valid_states:
        print(f"\nState {state_id+1}: {td.e[state_id]*27.211:.3f} eV")
        
        # 1. Transition density matrix
        if GENERATE_TRANSITION_DENSITY:
            dm_trans = calculate_transition_density_matrix(td, state_id)
            filename_trans = os.path.join(OUTPUT_DIR, f'transition_density_state{state_id+1}.cube')
            cubegen.density(mol, filename_trans, dm_trans, nx=nx, ny=ny, nz=nz)
            print(f"  ✓ Transition density: {filename_trans}")
        
        # 2. Excited state density
        if GENERATE_EXCITED_DENSITY:
            dm_excited = calculate_excited_state_density(td, state_id)
            filename_excited = os.path.join(OUTPUT_DIR, f'excited_state_density_state{state_id+1}.cube')
            cubegen.density(mol, filename_excited, dm_excited, nx=nx, ny=ny, nz=nz)
            print(f"  ✓ Excited state density: {filename_excited}")
        
        # 3. Density difference
        if GENERATE_DENSITY_DIFFERENCE:
            if not GENERATE_EXCITED_DENSITY:
                dm_excited = calculate_excited_state_density(td, state_id)
            
            # Get ground state density and convert CuPy to NumPy if needed
            dm_ground = mf.make_rdm1()
            if hasattr(dm_ground, 'get'):
                dm_ground = dm_ground.get()
            
            # For UKS, dm_ground is a tuple (alpha, beta), sum them
            if isinstance(dm_ground, tuple):
                dm_ground = dm_ground[0] + dm_ground[1]
            
            dm_diff = dm_excited - dm_ground
            filename_diff = os.path.join(OUTPUT_DIR, f'density_difference_state{state_id+1}.cube')
            cubegen.density(mol, filename_diff, dm_diff, nx=nx, ny=ny, nz=nz)
            print(f"  ✓ Density difference: {filename_diff}")
        
        # Quantitative verification for first state
        if state_id == 0 and GENERATE_HOMO_LUMO and GENERATE_TRANSITION_DENSITY:
            print("\n  " + "-"*66)
            print("  QUANTITATIVE VERIFICATION: S1 vs HOMO→LUMO")
            print("  " + "-"*66)
            
            # Get HOMO and LUMO indices (handle both RKS and UKS)
            if actual_spin == 1:
                mo_occ_ver = mf.mo_occ
                mo_coeff_ver = mf.mo_coeff
            else:
                mo_occ_ver = mf.mo_occ[0]
                mo_coeff_ver = mf.mo_coeff[0]
            
            # Convert CuPy to NumPy if needed
            if hasattr(mo_occ_ver, 'get'):
                mo_occ_ver = mo_occ_ver.get()
                mo_coeff_ver = mo_coeff_ver.get()
            
            homo_idx = np.where(mo_occ_ver > 0)[0][-1]
            lumo_idx = np.where(mo_occ_ver == 0)[0][0]
            
            # Calculate analytical HOMO→LUMO transition density
            homo_mo = mo_coeff_ver[:, homo_idx]
            lumo_mo = mo_coeff_ver[:, lumo_idx]
            t_homo_lumo = np.outer(homo_mo, lumo_mo) + np.outer(lumo_mo, homo_mo)
            
            # Get TDDFT transition density for S1
            dm_trans_s1 = calculate_transition_density_matrix(td, 0)
            
            # Calculate overlap/similarity (Frobenius inner product)
            overlap = np.sum(dm_trans_s1 * t_homo_lumo)
            norm_tddft = np.linalg.norm(dm_trans_s1)
            norm_homo_lumo = np.linalg.norm(t_homo_lumo)
            similarity = overlap / (norm_tddft * norm_homo_lumo)
            
            # Calculate HOMO→LUMO contribution from TDDFT amplitudes
            X, Y = td.xy[0]
            
            # For UKS, X and Y are tuples (Xa, Xb), (Ya, Yb)
            if actual_spin > 1:
                Xa, Xb = X
                Ya, Yb = Y
                # Convert CuPy to NumPy if needed
                if hasattr(Xa, 'get'):
                    Xa = Xa.get()
                    Ya = Ya.get()
                # Use only alpha part for verification
                # HOMO is last occupied (nocc_a-1), LUMO is first virtual (0)
                homo_lumo_amplitude = abs(Xa[-1, 0] + Ya[-1, 0])
                total_amplitude = np.linalg.norm(Xa + Ya)
            else:
                # Convert CuPy to NumPy if needed
                if hasattr(X, 'get'):
                    X = X.get()
                    Y = Y.get()
                nocc = X.shape[0]
                nvir = X.shape[1]
                # HOMO is index nocc-1 in occupied space, LUMO is index 0 in virtual space
                homo_lumo_amplitude = abs(X[nocc-1, 0] + Y[nocc-1, 0])
                total_amplitude = np.linalg.norm(X + Y)
            
            homo_lumo_weight = (homo_lumo_amplitude / total_amplitude)**2
            
            print(f"  Similarity (cosine): {similarity:.4f}")
            print(f"  HOMO→LUMO weight: {homo_lumo_weight:.4f} ({homo_lumo_weight*100:.1f}%)")
            print(f"  HOMO→LUMO amplitude: {homo_lumo_amplitude:.4f}")
            
            if similarity > 0.95 and homo_lumo_weight > 0.8:
                print("  ✓ S1 is STRONGLY dominated by HOMO→LUMO transition")
            elif similarity > 0.85 and homo_lumo_weight > 0.6:
                print("  ✓ S1 is MOSTLY a HOMO→LUMO transition")
            elif similarity > 0.70 and homo_lumo_weight > 0.4:
                print("  ⚠ S1 has SIGNIFICANT HOMO→LUMO character but mixed")
            else:
                print("  ⚠ S1 is NOT a pure HOMO→LUMO transition (multi-configurational)")
            
            print("  " + "-"*66)

print("="*70)

# ============================================================================
# 9A. GENERATE ORBITAL PAIR TRANSITION DENSITY CUBE FILES
# ============================================================================

if ENABLE_CONTRIBUTION_ANALYSIS and GENERATE_PAIR_CUBES and 'all_contributions' in locals():
    print("\n" + "="*70)
    print("GENERATING ORBITAL PAIR TRANSITION DENSITY CUBE FILES")
    print("="*70)
    
    for state_id in valid_contrib_states:
        if state_id not in all_contributions:
            continue
        
        contributions, _ = all_contributions[state_id]
        excitation_energy = td.e[state_id] * 27.211
        
        print(f"\nState {state_id+1} ({excitation_energy:.4f} eV):")
        
        pair_count = 0
        for rank, (occ_idx, vir_idx, weight, label) in enumerate(contributions, 1):
            if pair_count >= MAX_PAIRS_PER_STATE:
                break
            
            if weight < PAIR_CONTRIBUTION_CUTOFF:
                print(f"  Skipping {label} (contribution {weight*100:.2f}% < {PAIR_CONTRIBUTION_CUTOFF*100:.0f}%)")
                continue
            
            # Calculate transition density for this pair
            t_dm_pair = calculate_pair_transition_density(mf, occ_idx, vir_idx)
            
            # Debug: Check dimensions
            nao = mol.nao_nr()
            if t_dm_pair.shape[0] != nao or t_dm_pair.shape[1] != nao:
                print(f"  ⚠ Warning: Density matrix shape {t_dm_pair.shape} doesn't match AO basis {nao}x{nao}")
                print(f"  Skipping {label} - dimension mismatch")
                continue
            
            # Generate cube file
            labels, _ = get_orbital_labels(mf)
            occ_label = labels[occ_idx].replace('-', 'm').replace('+', 'p')
            vir_label = labels[vir_idx].replace('-', 'm').replace('+', 'p')
            
            filename = os.path.join(OUTPUT_DIR, 
                f'transition_pair_state{state_id+1}_{occ_label}_to_{vir_label}.cube')
            
            try:
                cubegen.density(mol, filename, t_dm_pair, nx=nx, ny=ny, nz=nz)
                print(f"  ✓ Rank {rank}: {label} ({weight*100:.2f}%) → {filename}")
                pair_count += 1
            except (AssertionError, ValueError) as e:
                print(f"  ✗ Failed to generate cube for {label}: {str(e)}")
                continue
    
    print("="*70)

# ============================================================================
# 10. CALCULATION SUMMARY
# ============================================================================

print("\n" + "="*70)
print("CALCULATION SUMMARY")
print("="*70)

print(f"\nMolecule: {'XYZ file: ' + XYZ_FILE if USE_XYZ else 'H2O test molecule'}")
print(f"Basis set: {BASIS_SET}")
print(f"Number of atoms: {mol.natm}")
print(f"Number of electrons: {mol.nelectron}")
print(f"Number of basis functions: {mol.nao}")
print(f"Molecular charge: {mol.charge}")
print(f"Spin multiplicity: {actual_spin} ({'closed-shell' if actual_spin == 1 else 'open-shell'})")

print(f"\nComputational settings:")
print(f"  DFT method: {dft_method}")
print(f"  XC functional: {XC_FUNCTIONAL}")
print(f"  Parallel threads: {NUM_THREADS if NUM_THREADS > 0 else 'auto-detected'}")
print(f"  TDDFT method: Full TDDFT (not TDA)")
print(f"  Number of excited states calculated: {NUM_EXCITED_STATES}")

print(f"\nGrid settings:")
if USE_GRID_RESOLUTION:
    print(f"  Mode: Fixed resolution")
    print(f"  Grid: {nx} × {ny} × {nz} = {nx*ny*nz:,} points")
else:
    print(f"  Mode: Box dimensions")
    print(f"  Margin: {BOX_MARGIN} Å, Spacing: {GRID_SPACING} Å")
    print(f"  Grid: {nx} × {ny} × {nz} = {nx*ny*nz:,} points")

print(f"\nOutput files generated:")
if GENERATE_HOMO_LUMO:
    print(f"  ✓ HOMO/LUMO orbitals (in {OUTPUT_DIR}/)")
if ENABLE_NTO_ANALYSIS and valid_nto_states:
    print(f"  ✓ NTO analysis for {len(valid_nto_states)} state(s)")
if valid_states:
    print(f"  ✓ Cube files for {len(valid_states)} state(s):")
    if GENERATE_TRANSITION_DENSITY:
        print(f"    - Transition density matrices")
    if GENERATE_EXCITED_DENSITY:
        print(f"    - Excited state densities")
    if GENERATE_DENSITY_DIFFERENCE:
        print(f"    - Density differences")

print(f"\nAll output files saved to: {OUTPUT_DIR}/")

# ============================================================================
# VISUALIZATION GUIDE (Commented out to reduce log clutter)
# ============================================================================
# print("\n" + "="*70)
# print("VISUALIZATION GUIDE")
# print("="*70)
# print("""
# Three types of cube files generated:
# 
# 1. transition_density_state*.cube
#    - Represents the electronic transition between ground and excited state
#    - Used for calculating transition dipole moments
#    - Visualize with isovalues ±0.002
# 
# 2. excited_state_density_state*.cube
#    - Total electron density in the excited state
#    - Compare with ground state density
# 
# 3. density_difference_state*.cube
#    - Change in electron density (excited - ground)
#    - Red/positive: electron accumulation
#    - Blue/negative: electron depletion
#    - Recommended for visualization
# 
# VMD visualization:
#   vmd output/density_difference_state1.cube
#   Graphics > Representations > Drawing Method: Isosurface
#   - Rep 1: Isovalue = +0.002, Color = Red (electron gain)
#   - Rep 2: Isovalue = -0.002, Color = Blue (electron loss)
# 
# Jmol visualization:
#   isosurface ID "surf1" cutoff  0.002 output/density_difference_state1.cube
#   isosurface ID "surf2" cutoff -0.002 output/density_difference_state1.cube
# 
# HOMO/LUMO verification:
#   Compare transition_density_state1.cube with HOMO.cube and LUMO.cube
#   to verify that S1 corresponds to a HOMO→LUMO transition.
# 
# NTO visualization:
#   Open output/nto_state_*.molden files in Jmol, Avogadro, or VMD
# """)
# print("="*70)

print("\n✓ Calculation completed successfully!")
print(f"✓ All files saved to: {OUTPUT_DIR}/\n")
