#!/usr/bin/env python
'''
Analyze Orbital Contributions to TDDFT Transitions

This script analyzes which orbital pairs contribute to each excited state
and generates transition density maps for individual orbital pair contributions.

Similar to GPAW's transition contribution analysis:
https://gpaw.readthedocs.io/tutorialsexercises/opticalresponse/tddft/lcaotddft.html

Features:
- Analyze orbital pair contributions (i→a) for each excited state
- Show percentage contribution of each pair
- Generate transition density cube files for dominant pairs
- Create detailed contribution tables
'''

import numpy as np
import os
from pyscf import gto, dft, tddft, lib
from pyscf.tools import cubegen, molden

# ============================================================================
# CONFIGURATION
# ============================================================================

# --- Molecule Selection ---
USE_XYZ = True
XYZ_FILE = 'PTCDA.xyz'
BASIS_SET = '6-31g'

# --- DFT/TDDFT Settings ---
XC_FUNCTIONAL = 'b3lyp'
NUM_EXCITED_STATES = 10

# --- Analysis Settings ---
STATES_TO_ANALYZE = [0, 1, 2]  # Which states to analyze (0-indexed)
CONTRIBUTION_THRESHOLD = 0.01  # Show contributions > 1%
TOP_N_CONTRIBUTIONS = 10       # Show top N orbital pairs

# --- Cube File Generation ---
GENERATE_PAIR_CUBES = True     # Generate cube files for orbital pairs
MAX_PAIRS_PER_STATE = 3        # Generate cubes for top N pairs per state
PAIR_CONTRIBUTION_CUTOFF = 0.05  # Only generate cubes for pairs > 5%

# --- Grid Settings ---
USE_GRID_RESOLUTION = False
GRID_RESOLUTION = [80, 80, 80]
BOX_MARGIN = 4.0  # Angstrom
GRID_SPACING = 0.2  # Angstrom

# --- Output ---
OUTPUT_DIR = 'transition_analysis'
SAVE_CONTRIBUTION_TABLE = True

# --- Parallel Settings ---
NUM_THREADS = 0  # 0 = auto-detect

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def setup_molecule():
    """Setup molecule from XYZ file or test molecule."""
    mol = gto.Mole()
    
    if USE_XYZ:
        mol.atom = XYZ_FILE
    else:
        mol.atom = 'H 0 0 0; F 0 0 1.1'
    
    mol.basis = BASIS_SET
    mol.build()
    return mol

def calculate_grid_parameters(mol, use_resolution=True, resolution=None, 
                              box_margin=4.0, grid_spacing=0.2):
    """Calculate grid parameters for cube files."""
    coords = mol.atom_coords()
    mol_min = coords.min(axis=0)
    mol_max = coords.max(axis=0)
    mol_size = mol_max - mol_min
    
    box_info = {'mol_size': mol_size}
    
    if use_resolution:
        nx, ny, nz = resolution
        box_info['grid_resolution'] = resolution
    else:
        box_size = mol_size + 2 * box_margin
        nx = int(np.ceil(box_size[0] / grid_spacing))
        ny = int(np.ceil(box_size[1] / grid_spacing))
        nz = int(np.ceil(box_size[2] / grid_spacing))
        box_info['grid_spacing'] = [grid_spacing] * 3
        box_info['box_size'] = box_size
    
    box_info['total_points'] = nx * ny * nz
    return nx, ny, nz, box_info

def get_orbital_labels(mf):
    """Get orbital labels (HOMO-n, LUMO+n)."""
    mo_occ = mf.mo_occ
    homo_idx = np.where(mo_occ > 0)[0][-1]
    
    labels = []
    for i in range(len(mo_occ)):
        if i <= homo_idx:
            # Occupied orbital
            offset = homo_idx - i
            if offset == 0:
                labels.append('HOMO')
            else:
                labels.append(f'HOMO-{offset}')
        else:
            # Virtual orbital
            offset = i - homo_idx - 1
            if offset == 0:
                labels.append('LUMO')
            else:
                labels.append(f'LUMO+{offset}')
    
    return labels, homo_idx

def analyze_transition_contributions(td, state_id, mf, threshold=0.01, top_n=10):
    """
    Analyze orbital pair contributions to a specific excited state.
    
    Returns:
        contributions: list of (occ_idx, vir_idx, weight, label) sorted by weight
    """
    X, Y = td.xy[state_id]
    nocc, nvir = X.shape
    
    # For full TDDFT, the contribution is from (X + Y)
    # For TDA, Y = 0
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
    contributions = [(occ, vir, w/total_weight, label) 
                     for occ, vir, w, label in contributions]
    
    # Return top N
    return contributions[:top_n], total_weight

def calculate_pair_transition_density(mf, occ_idx, vir_idx):
    """
    Calculate transition density matrix for a single orbital pair i→a.
    
    T_μν = C_μ^i × C_ν^a + C_μ^a × C_ν^i
    """
    mo_coeff = mf.mo_coeff
    
    occ_mo = mo_coeff[:, occ_idx]
    vir_mo = mo_coeff[:, vir_idx]
    
    # Transition density matrix for this pair
    t_dm = np.outer(occ_mo, vir_mo) + np.outer(vir_mo, occ_mo)
    
    return t_dm

def print_contribution_table(state_id, excitation_energy, contributions, total_weight):
    """Print a formatted table of orbital contributions."""
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
    print(f"{'='*70}")

def save_contribution_table_to_file(filename, all_contributions, td):
    """Save contribution tables to a text file."""
    with open(filename, 'w') as f:
        f.write("="*70 + "\n")
        f.write("ORBITAL PAIR CONTRIBUTIONS TO EXCITED STATES\n")
        f.write("="*70 + "\n\n")
        
        for state_id, (contributions, total_weight) in all_contributions.items():
            excitation_energy = td.e[state_id] * 27.211  # Convert to eV
            
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

# ============================================================================
# MAIN SCRIPT
# ============================================================================

def main():
    print("="*70)
    print("TRANSITION CONTRIBUTION ANALYSIS")
    print("="*70)
    
    # Setup parallel computation
    if NUM_THREADS > 0:
        lib.num_threads(NUM_THREADS)
        print(f"\nUsing {NUM_THREADS} threads for parallel computation")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_DIR}/")
    
    # Setup molecule
    print("\n" + "="*70)
    print("MOLECULE SETUP")
    print("="*70)
    mol = setup_molecule()
    print(f"Molecule: {XYZ_FILE if USE_XYZ else 'H-F test molecule'}")
    print(f"Basis set: {BASIS_SET}")
    print(f"Number of atoms: {mol.natm}")
    print(f"Number of electrons: {mol.nelectron}")
    print(f"Number of basis functions: {mol.nao}")
    
    # Ground state DFT
    print("\n" + "="*70)
    print("GROUND STATE DFT")
    print("="*70)
    print(f"XC functional: {XC_FUNCTIONAL}")
    
    mf = dft.RKS(mol)
    mf.xc = XC_FUNCTIONAL
    mf.kernel()
    
    print(f"Ground state energy: {mf.e_tot:.6f} a.u.")
    
    # Get orbital information
    mo_occ = mf.mo_occ
    homo_idx = np.where(mo_occ > 0)[0][-1]
    lumo_idx = np.where(mo_occ == 0)[0][0]
    nocc = np.sum(mo_occ > 0)
    nvir = len(mo_occ) - nocc
    
    print(f"Number of occupied orbitals: {nocc}")
    print(f"Number of virtual orbitals: {nvir}")
    print(f"HOMO index: {homo_idx}")
    print(f"LUMO index: {lumo_idx}")
    
    # TDDFT calculation
    print("\n" + "="*70)
    print("TDDFT CALCULATION")
    print("="*70)
    
    td = tddft.TDDFT(mf)
    td.nstates = NUM_EXCITED_STATES
    print(f"Calculating {NUM_EXCITED_STATES} excited states...")
    td.kernel()
    
    # Calculate grid parameters
    print("\n" + "="*70)
    print("GRID PARAMETERS")
    print("="*70)
    
    if USE_GRID_RESOLUTION:
        nx, ny, nz, box_info = calculate_grid_parameters(
            mol, use_resolution=True, resolution=GRID_RESOLUTION
        )
        print(f"Using fixed grid resolution: {nx} × {ny} × {nz}")
    else:
        nx, ny, nz, box_info = calculate_grid_parameters(
            mol, use_resolution=False, box_margin=BOX_MARGIN, 
            grid_spacing=GRID_SPACING
        )
        print(f"Using box dimensions with margin: {BOX_MARGIN} Å")
        print(f"Grid spacing: {GRID_SPACING} Å")
        print(f"Calculated grid resolution: {nx} × {ny} × {nz}")
    
    print(f"Total grid points: {box_info['total_points']:,}")
    
    # Analyze contributions for selected states
    print("\n" + "="*70)
    print("ANALYZING ORBITAL CONTRIBUTIONS")
    print("="*70)
    
    all_contributions = {}
    
    for state_id in STATES_TO_ANALYZE:
        if state_id >= td.nstates:
            print(f"\nWarning: State {state_id+1} not available (only {td.nstates} states calculated)")
            continue
        
        excitation_energy = td.e[state_id] * 27.211  # Convert to eV
        
        contributions, total_weight = analyze_transition_contributions(
            td, state_id, mf, 
            threshold=CONTRIBUTION_THRESHOLD,
            top_n=TOP_N_CONTRIBUTIONS
        )
        
        all_contributions[state_id] = (contributions, total_weight)
        
        # Print contribution table
        print_contribution_table(state_id, excitation_energy, contributions, total_weight)
    
    # Save contribution tables to file
    if SAVE_CONTRIBUTION_TABLE:
        table_file = os.path.join(OUTPUT_DIR, 'contribution_tables.txt')
        save_contribution_table_to_file(table_file, all_contributions, td)
        print(f"\n✓ Contribution tables saved to: {table_file}")
    
    # Generate cube files for dominant orbital pairs
    if GENERATE_PAIR_CUBES:
        print("\n" + "="*70)
        print("GENERATING ORBITAL PAIR TRANSITION DENSITY CUBE FILES")
        print("="*70)
        
        for state_id in STATES_TO_ANALYZE:
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
                
                # Generate cube file
                labels, _ = get_orbital_labels(mf)
                occ_label = labels[occ_idx].replace('-', 'm').replace('+', 'p')
                vir_label = labels[vir_idx].replace('-', 'm').replace('+', 'p')
                
                filename = os.path.join(OUTPUT_DIR, 
                    f'transition_pair_state{state_id+1}_{occ_label}_to_{vir_label}.cube')
                
                cubegen.density(mol, filename, t_dm_pair, nx=nx, ny=ny, nz=nz)
                
                print(f"  ✓ Rank {rank}: {label} ({weight*100:.2f}%) → {filename}")
                pair_count += 1
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    print(f"\nAnalyzed {len(all_contributions)} excited states")
    print(f"Output directory: {OUTPUT_DIR}/")
    
    print("\nGenerated files:")
    if SAVE_CONTRIBUTION_TABLE:
        print(f"  1. contribution_tables.txt - Detailed contribution tables")
    
    if GENERATE_PAIR_CUBES:
        total_cubes = sum(
            min(MAX_PAIRS_PER_STATE, 
                sum(1 for _, _, w, _ in all_contributions[sid][0] 
                    if w >= PAIR_CONTRIBUTION_CUTOFF))
            for sid in all_contributions
        )
        print(f"  2. {total_cubes} orbital pair transition density cube files")
    
    print("\nVisualization:")
    print("  Load cube files in VMD/Jmol to see individual orbital pair contributions")
    print("  Compare with full TDDFT transition density to verify")
    
    print("\nKey findings:")
    for state_id in STATES_TO_ANALYZE:
        if state_id not in all_contributions:
            continue
        
        contributions, _ = all_contributions[state_id]
        excitation_energy = td.e[state_id] * 27.211
        
        if contributions:
            top_contrib = contributions[0]
            print(f"  State {state_id+1} ({excitation_energy:.4f} eV): "
                  f"{top_contrib[3]} dominates with {top_contrib[2]*100:.1f}%")
    
    print("\n" + "="*70)
    print("DONE")
    print("="*70)

if __name__ == "__main__":
    main()
