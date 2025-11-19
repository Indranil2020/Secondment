#!/usr/bin/env python
'''
Accurate Transition Density Matrix and Transition Dipole Moment Calculation
Based on official PySCF examples:
- examples/tddft/22-density.py (density matrices)
- examples/1-advanced/030-transition_dipole.py (transition dipoles)
- examples/tddft/01-nto_analysis.py (NTO analysis)
'''

from pyscf import gto, dft, tddft
from pyscf.tools import cubegen, molden
import numpy as np
from functools import reduce

# ============================================================================
# 1. MOLECULE DEFINITION
# ============================================================================

# Option A: Simple molecule (H2O) for testing
def create_h2o_molecule():
    mol = gto.M(
        atom = '''
        O  0.0000  0.0000  0.1173
        H  0.0000  0.7572 -0.4692
        H  0.0000 -0.7572 -0.4692
        ''',
        basis = '6-31g',
        verbose = 4
    )
    return mol

# Option B: Load from XYZ file (PTCDA)
def create_molecule_from_xyz(xyz_file='PTCDA.xyz', basis='6-31g'):
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

# Choose molecule
# USE_XYZ = False  # Set to True to use PTCDA.xyz
USE_XYZ = True  # Set to True to use PTCDA.xyz
if USE_XYZ:
    mol = create_molecule_from_xyz('PTCDA.xyz', basis='6-31g')
else:
    mol = create_h2o_molecule()

# ============================================================================
# 2. GROUND STATE DFT CALCULATION
# ============================================================================

print("\n" + "="*70)
print("GROUND STATE DFT CALCULATION")
print("="*70)

mf = dft.RKS(mol)
mf.xc = 'b3lyp'  # B3LYP functional
mf.kernel()

print(f"Ground state energy: {mf.e_tot:.6f} a.u.")

# ============================================================================
# 3. TDDFT CALCULATION
# ============================================================================

print("\n" + "="*70)
print("TDDFT CALCULATION")
print("="*70)

# Use full TDDFT (not TDA) for more accurate transition dipoles
td = tddft.TDDFT(mf)
td.nstates = 10  # Calculate first 10 excited states
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

print("\n" + "="*70)
print("NATURAL TRANSITION ORBITAL ANALYSIS")
print("="*70)

# Perform NTO analysis for first few states
nto_states = min(3, td.nstates)  # Analyze first 3 states

for i in range(nto_states):
    print(f"\nState {i+1} ({td.e[i]*27.211:.3f} eV):")
    weights, nto_coeff = td.get_nto(state=i+1, verbose=4)
    
    # Save NTO orbitals to molden format for visualization
    molden_file = f'nto_state_{i+1}.molden'
    molden.from_mo(mol, molden_file, nto_coeff)
    print(f"  NTO orbitals saved to: {molden_file}")

print("="*70)

# ============================================================================
# 7. GENERATE CUBE FILES FOR VISUALIZATION
# ============================================================================

print("\n" + "="*70)
print("GENERATING CUBE FILES")
print("="*70)

# Select states to visualize
states_to_visualize = [0, 1, 2]  # First three states (0-indexed)
states_to_visualize = [s for s in states_to_visualize if s < td.nstates]

# Grid resolution
nx, ny, nz = 80, 80, 80

for state_id in states_to_visualize:
    print(f"\nState {state_id+1}: {td.e[state_id]*27.211:.3f} eV")
    
    # 1. Transition density matrix (most important for transition moments)
    dm_trans = calculate_transition_density_matrix(td, state_id)
    filename_trans = f'transition_density_state{state_id+1}.cube'
    cubegen.density(mol, filename_trans, dm_trans, nx=nx, ny=ny, nz=nz)
    print(f"  ✓ Transition density: {filename_trans}")
    
    # 2. Excited state density
    dm_excited = calculate_excited_state_density(td, state_id)
    filename_excited = f'excited_state_density_state{state_id+1}.cube'
    cubegen.density(mol, filename_excited, dm_excited, nx=nx, ny=ny, nz=nz)
    print(f"  ✓ Excited state density: {filename_excited}")
    
    # 3. Density difference (excited - ground)
    dm_diff = dm_excited - mf.make_rdm1()
    filename_diff = f'density_difference_state{state_id+1}.cube'
    cubegen.density(mol, filename_diff, dm_diff, nx=nx, ny=ny, nz=nz)
    print(f"  ✓ Density difference: {filename_diff}")

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
  vmd density_difference_state1.cube
  Graphics > Representations > Drawing Method: Isosurface
  - Rep 1: Isovalue = +0.002, Color = Red (electron gain)
  - Rep 2: Isovalue = -0.002, Color = Blue (electron loss)

Jmol visualization:
  isosurface ID "surf1" cutoff  0.002 density_difference_state1.cube
  isosurface ID "surf2" cutoff -0.002 density_difference_state1.cube

NTO visualization:
  Open nto_state_*.molden files in Jmol, Avogadro, or VMD
""")
print("="*70)

print("\n✓ Calculation completed successfully!")
