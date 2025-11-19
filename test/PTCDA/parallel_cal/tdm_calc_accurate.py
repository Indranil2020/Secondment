#!/usr/bin/env python
'''
Accurate Transition Density Matrix and Transition Dipole Moment Calculation
Based on official PySCF examples:
- examples/tddft/22-density.py (density matrices)
- examples/1-advanced/030-transition_dipole.py (transition dipoles)
- examples/tddft/01-nto_analysis.py (NTO analysis)

Features:
- Parallel calculation support
- Configurable grid size and box dimensions
- Selective state output
- HOMO/LUMO cube file generation
'''

from pyscf import gto, dft, tddft, lib
from pyscf.tools import cubegen, molden
import numpy as np
from functools import reduce
import os

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE SETTINGS
# ============================================================================

# --- Parallel Calculation Settings ---
ENABLE_PARALLEL = True  # Enable/disable parallel computation
NUM_THREADS = 0  # Number of CPU threads (set to 0 for auto-detect)

# --- Molecule Selection ---
USE_XYZ = True  # True: use XYZ file, False: use H2O test molecule
XYZ_FILE = 'PTCDA.xyz'  # Path to XYZ file
BASIS_SET = '6-31g'  # Basis set

# --- DFT/TDDFT Settings ---
# Note: TDDFT uses the same basis set and XC functional as ground state DFT
XC_FUNCTIONAL = 'b3lyp'  # Exchange-correlation functional
# Common options:
#   'b3lyp'    - B3LYP (hybrid, good general purpose)
#   'pbe0'     - PBE0 (hybrid, good for excited states)
#   'cam-b3lyp' - CAM-B3LYP (range-separated, good for charge transfer)
#   'wb97x-d'  - ωB97X-D (range-separated with dispersion)
#   'pbe'      - PBE (GGA, faster but less accurate)
#   'blyp'     - BLYP (GGA)

NUM_EXCITED_STATES = 10  # Total number of excited states to calculate

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
USE_GRID_RESOLUTION = True
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
OUTPUT_DIR = 'output'  # Directory for output files (created if doesn't exist)

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
        verbose = 4
    )
    return mol

def create_molecule_from_xyz(xyz_file, basis):
    """
    Load molecule from XYZ file.
    Note: XYZ file should be in standard format (first line: number of atoms)
    """
    mol = gto.M(
        atom = xyz_file,
        basis = basis,
        verbose = 4
    )
    return mol

print("\n" + "="*70)
print("MOLECULE SETUP")
print("="*70)

if USE_XYZ:
    print(f"Loading molecule from: {XYZ_FILE}")
    print(f"Basis set: {BASIS_SET}")
    mol = create_molecule_from_xyz(XYZ_FILE, BASIS_SET)
else:
    print("Using H2O test molecule")
    print(f"Basis set: {BASIS_SET}")
    mol = create_h2o_molecule()

print(f"Number of atoms: {mol.natm}")
print(f"Number of electrons: {mol.nelectron}")
print(f"Number of basis functions: {mol.nao}")

# ============================================================================
# 2. GROUND STATE DFT CALCULATION
# ============================================================================

print("\n" + "="*70)
print("GROUND STATE DFT CALCULATION")
print("="*70)
print(f"XC functional: {XC_FUNCTIONAL}")
print(f"Basis set: {BASIS_SET}")

mf = dft.RKS(mol)
mf.xc = XC_FUNCTIONAL
mf.kernel()

print(f"Ground state energy: {mf.e_tot:.6f} a.u.")

# ============================================================================
# 3. TDDFT CALCULATION
# ============================================================================

print("\n" + "="*70)
print("TDDFT CALCULATION")
print("="*70)
print(f"Note: TDDFT inherits XC functional ({XC_FUNCTIONAL}) and basis set ({BASIS_SET}) from ground state")

# Use full TDDFT (not TDA) for more accurate transition dipoles
td = tddft.TDDFT(mf)
td.nstates = NUM_EXCITED_STATES
print(f"Calculating {NUM_EXCITED_STATES} excited states...")
td.kernel()
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
    Calculate transition dipole moment between ground state and excited state.
    
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
    
    # Get MO coefficients
    mo_coeff = td._scf.mo_coeff
    mo_occ = td._scf.mo_occ
    
    # Occupied and virtual orbitals
    orbo = mo_coeff[:, mo_occ > 0]
    orbv = mo_coeff[:, mo_occ == 0]
    
    # Construct transition density matrix in AO basis
    # For RHF/RKS: T = C_occ * (X + Y) * C_vir^T
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
    
    orbo = mo_coeff[:, mo_occ > 0]
    orbv = mo_coeff[:, mo_occ == 0]
    
    nocc = orbo.shape[1]
    nvir = orbv.shape[1]
    
    # Transition density in MO basis
    t_dm1_mo = np.zeros((mo_coeff.shape[1], mo_coeff.shape[1]))
    t_dm1_mo[:nocc, nocc:] = (X + Y).reshape(nocc, nvir)
    t_dm1_mo[nocc:, :nocc] = (X + Y).reshape(nocc, nvir).T
    
    # Transform to AO basis
    t_dm1_ao = reduce(np.dot, (mo_coeff, t_dm1_mo, mo_coeff.T))
    
    return t_dm1_ao

def calculate_excited_state_density(td, state_id):
    """
    Calculate excited state density matrix.
    Based on PySCF examples/tddft/22-density.py
    
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
    
    # For TDDFT, need to consider both X and Y
    # Density matrix changes in MO basis
    nocc = X.shape[0]
    
    # Occupied-occupied and virtual-virtual blocks
    dm_oo = -np.einsum('ia,ka->ik', X.conj(), X)
    dm_oo -= np.einsum('ia,ka->ik', Y.conj(), Y)
    
    dm_vv = np.einsum('ia,ic->ac', X, X.conj())
    dm_vv += np.einsum('ia,ic->ac', Y, Y.conj())
    
    # Start with ground state density in MO basis
    mf = td._scf
    dm = np.diag(mf.mo_occ)
    
    # Add TDDFT contribution
    dm[:nocc, :nocc] += dm_oo * 2
    dm[nocc:, nocc:] += dm_vv * 2
    
    # Transform to AO basis
    mo = mf.mo_coeff
    dm_excited = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    
    return dm_excited

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
            weights, nto_coeff = td.get_nto(state=i+1, verbose=4)
            
            # Save NTO orbitals to molden format for visualization
            molden_file = os.path.join(OUTPUT_DIR, f'nto_state_{i+1}.molden')
            molden.from_mo(mol, molden_file, nto_coeff)
            print(f"  NTO orbitals saved to: {molden_file}")
    
    print("="*70)
else:
    print("\nNTO analysis disabled.")

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
    
    # Find HOMO and LUMO indices
    mo_occ = mf.mo_occ
    homo_idx = np.where(mo_occ > 0)[0][-1]  # Last occupied orbital
    lumo_idx = np.where(mo_occ == 0)[0][0]  # First unoccupied orbital
    
    print(f"HOMO index: {homo_idx} (MO {homo_idx+1})")
    print(f"LUMO index: {lumo_idx} (MO {lumo_idx+1})")
    print(f"HOMO energy: {mf.mo_energy[homo_idx]*27.211:.3f} eV")
    print(f"LUMO energy: {mf.mo_energy[lumo_idx]*27.211:.3f} eV")
    print(f"HOMO-LUMO gap: {(mf.mo_energy[lumo_idx] - mf.mo_energy[homo_idx])*27.211:.3f} eV")
    
    # Generate HOMO cube file
    homo_file = os.path.join(OUTPUT_DIR, 'HOMO.cube')
    cubegen.orbital(mol, homo_file, mf.mo_coeff[:, homo_idx], nx=nx, ny=ny, nz=nz)
    print(f"\n  ✓ HOMO orbital: {homo_file}")
    
    # Generate LUMO cube file
    lumo_file = os.path.join(OUTPUT_DIR, 'LUMO.cube')
    cubegen.orbital(mol, lumo_file, mf.mo_coeff[:, lumo_idx], nx=nx, ny=ny, nz=nz)
    print(f"  ✓ LUMO orbital: {lumo_file}")
    
    # Generate HOMO-1 and LUMO+1 for additional verification
    if homo_idx > 0:
        homo1_file = os.path.join(OUTPUT_DIR, 'HOMO-1.cube')
        cubegen.orbital(mol, homo1_file, mf.mo_coeff[:, homo_idx-1], nx=nx, ny=ny, nz=nz)
        print(f"  ✓ HOMO-1 orbital: {homo1_file}")
    
    if lumo_idx < len(mo_occ) - 1:
        lumo1_file = os.path.join(OUTPUT_DIR, 'LUMO+1.cube')
        cubegen.orbital(mol, lumo1_file, mf.mo_coeff[:, lumo_idx+1], nx=nx, ny=ny, nz=nz)
        print(f"  ✓ LUMO+1 orbital: {lumo1_file}")
    
    print("\nVerification tip:")
    print("For the first excited state (S1), check if the transition density")
    print("resembles a HOMO→LUMO transition by comparing:")
    print("  - transition_density_state1.cube")
    print("  - HOMO.cube (electron depletion)")
    print("  - LUMO.cube (electron accumulation)")
    
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
            dm_diff = dm_excited - mf.make_rdm1()
            filename_diff = os.path.join(OUTPUT_DIR, f'density_difference_state{state_id+1}.cube')
            cubegen.density(mol, filename_diff, dm_diff, nx=nx, ny=ny, nz=nz)
            print(f"  ✓ Density difference: {filename_diff}")

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

print(f"\nComputational settings:")
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

print("\n" + "="*70)
print("VISUALIZATION GUIDE")
print("="*70)
print("""
Three types of cube files generated:

1. transition_density_state*.cube
   - Represents the electronic transition between ground and excited state
   - Used for calculating transition dipole moments
   - Visualize with isovalues ±0.002

2. excited_state_density_state*.cube
   - Total electron density in the excited state
   - Compare with ground state density

3. density_difference_state*.cube
   - Change in electron density (excited - ground)
   - Red/positive: electron accumulation
   - Blue/negative: electron depletion
   - Recommended for visualization

VMD visualization:
  vmd output/density_difference_state1.cube
  Graphics > Representations > Drawing Method: Isosurface
  - Rep 1: Isovalue = +0.002, Color = Red (electron gain)
  - Rep 2: Isovalue = -0.002, Color = Blue (electron loss)

Jmol visualization:
  isosurface ID "surf1" cutoff  0.002 output/density_difference_state1.cube
  isosurface ID "surf2" cutoff -0.002 output/density_difference_state1.cube

HOMO/LUMO verification:
  Compare transition_density_state1.cube with HOMO.cube and LUMO.cube
  to verify that S1 corresponds to a HOMO→LUMO transition.

NTO visualization:
  Open output/nto_state_*.molden files in Jmol, Avogadro, or VMD
""")
print("="*70)

print("\n✓ Calculation completed successfully!")
print(f"✓ All files saved to: {OUTPUT_DIR}/\n")
