#!/usr/bin/env python
'''
Compare HOMO×LUMO Multiplication vs TDDFT Transition Density

This script demonstrates the difference between:
1. Approximate method: ρ_approx(r) = φ_HOMO(r) × φ_LUMO(r)
2. Exact method: ρ_exact(r) from TDDFT transition density matrix

Scientific Note:
The approximation ρ(r) ≈ φ_HOMO(r) × φ_LUMO(r) is qualitatively correct
but quantitatively different from the exact TDDFT transition density.
The exact method properly accounts for basis function overlap and is
computed as: ρ(r) = Σ_μν T_μν χ_μ(r) χ_ν(r)
'''

import numpy as np
import os
import sys

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directory containing cube files (from tdm_calc_accurate.py output)
CUBE_DIR = 'output'

# Input cube files
HOMO_CUBE = os.path.join(CUBE_DIR, 'HOMO.cube')
LUMO_CUBE = os.path.join(CUBE_DIR, 'LUMO.cube')
TDDFT_TRANSITION_CUBE = os.path.join(CUBE_DIR, 'transition_density_state1.cube')
ANALYTICAL_CUBE = os.path.join(CUBE_DIR, 'transition_HOMO_LUMO_analytical.cube')

# Output cube file
OUTPUT_CUBE = os.path.join(CUBE_DIR, 'transition_HOMO_LUMO_multiplied.cube')

# Comparison settings
GENERATE_DIFFERENCE_CUBES = True  # Generate difference cube files
CALCULATE_STATISTICS = True       # Calculate detailed statistics

# ============================================================================
# CUBE FILE I/O FUNCTIONS
# ============================================================================

def read_cube_file(filename):
    """
    Read a Gaussian cube file.
    
    Returns:
        atoms: list of (atomic_number, charge, x, y, z)
        origin: (x, y, z) origin of the grid
        grid_vectors: [(nx, dx_x, dx_y, dx_z), (ny, dy_x, dy_y, dy_z), (nz, dz_x, dz_y, dz_z)]
        data: 3D numpy array of density values
    """
    with open(filename, 'r') as f:
        # Skip first two comment lines
        comment1 = f.readline()
        comment2 = f.readline()
        
        # Read number of atoms and origin
        line = f.readline().split()
        natoms = int(line[0])
        origin = np.array([float(line[1]), float(line[2]), float(line[3])])
        
        # Read grid dimensions and vectors
        grid_vectors = []
        for i in range(3):
            line = f.readline().split()
            n = int(line[0])
            vec = np.array([float(line[1]), float(line[2]), float(line[3])])
            grid_vectors.append((n, vec))
        
        # Read atom information
        atoms = []
        for i in range(abs(natoms)):
            line = f.readline().split()
            atomic_num = int(line[0])
            charge = float(line[1])
            x, y, z = float(line[2]), float(line[3]), float(line[4])
            atoms.append((atomic_num, charge, x, y, z))
        
        # Read volumetric data
        nx, ny, nz = grid_vectors[0][0], grid_vectors[1][0], grid_vectors[2][0]
        data = np.zeros((nx, ny, nz))
        
        values = []
        for line in f:
            values.extend([float(x) for x in line.split()])
        
        # Reshape data
        data = np.array(values).reshape((nx, ny, nz))
        
    return atoms, origin, grid_vectors, data

def write_cube_file(filename, atoms, origin, grid_vectors, data, comment1="", comment2=""):
    """
    Write a Gaussian cube file.
    """
    nx, ny, nz = data.shape
    
    with open(filename, 'w') as f:
        # Write comments
        f.write(f"{comment1}\n")
        f.write(f"{comment2}\n")
        
        # Write number of atoms and origin
        f.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        
        # Write grid vectors
        for n, vec in grid_vectors:
            f.write(f"{n:5d} {vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")
        
        # Write atoms
        for atom in atoms:
            atomic_num, charge, x, y, z = atom
            f.write(f"{atomic_num:5d} {charge:12.6f} {x:12.6f} {y:12.6f} {z:12.6f}\n")
        
        # Write volumetric data (6 values per line)
        values = data.flatten()
        for i in range(0, len(values), 6):
            line_values = values[i:i+6]
            f.write(" ".join([f"{v:13.5E}" for v in line_values]) + "\n")

# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================

def calculate_similarity(data1, data2):
    """
    Calculate cosine similarity between two 3D arrays.
    
    similarity = <data1, data2> / (||data1|| * ||data2||)
    """
    dot_product = np.sum(data1 * data2)
    norm1 = np.linalg.norm(data1)
    norm2 = np.linalg.norm(data2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    return dot_product / (norm1 * norm2)

def calculate_statistics(data1, data2, name1="Data1", name2="Data2"):
    """
    Calculate detailed statistics comparing two datasets.
    """
    print(f"\n{'='*70}")
    print(f"STATISTICAL COMPARISON: {name1} vs {name2}")
    print(f"{'='*70}")
    
    # Basic statistics
    print(f"\n{name1}:")
    print(f"  Min:  {np.min(data1):12.6e}")
    print(f"  Max:  {np.max(data1):12.6e}")
    print(f"  Mean: {np.mean(data1):12.6e}")
    print(f"  Std:  {np.std(data1):12.6e}")
    print(f"  Norm: {np.linalg.norm(data1):12.6e}")
    
    print(f"\n{name2}:")
    print(f"  Min:  {np.min(data2):12.6e}")
    print(f"  Max:  {np.max(data2):12.6e}")
    print(f"  Mean: {np.mean(data2):12.6e}")
    print(f"  Std:  {np.std(data2):12.6e}")
    print(f"  Norm: {np.linalg.norm(data2):12.6e}")
    
    # Difference statistics
    diff = data1 - data2
    print(f"\nDifference ({name1} - {name2}):")
    print(f"  Min:  {np.min(diff):12.6e}")
    print(f"  Max:  {np.max(diff):12.6e}")
    print(f"  Mean: {np.mean(diff):12.6e}")
    print(f"  Std:  {np.std(diff):12.6e}")
    print(f"  RMS:  {np.sqrt(np.mean(diff**2)):12.6e}")
    
    # Relative difference
    # Avoid division by zero
    mask = np.abs(data2) > 1e-10
    if np.any(mask):
        rel_diff = np.abs(diff[mask] / data2[mask])
        print(f"\nRelative Difference (where |{name2}| > 1e-10):")
        print(f"  Mean: {np.mean(rel_diff):12.6e}")
        print(f"  Max:  {np.max(rel_diff):12.6e}")
    
    # Similarity metrics
    similarity = calculate_similarity(data1, data2)
    print(f"\nSimilarity Metrics:")
    print(f"  Cosine similarity: {similarity:8.6f}")
    
    # Correlation coefficient
    corr = np.corrcoef(data1.flatten(), data2.flatten())[0, 1]
    print(f"  Correlation coef:  {corr:8.6f}")
    
    # R² (coefficient of determination)
    ss_res = np.sum((data1 - data2)**2)
    ss_tot = np.sum((data1 - np.mean(data1))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
    print(f"  R² score:          {r_squared:8.6f}")
    
    print(f"\n{'='*70}")
    
    return {
        'similarity': similarity,
        'correlation': corr,
        'r_squared': r_squared,
        'rms_diff': np.sqrt(np.mean(diff**2))
    }

# ============================================================================
# MAIN SCRIPT
# ============================================================================

def main():
    print("="*70)
    print("HOMO×LUMO MULTIPLICATION vs TDDFT TRANSITION DENSITY")
    print("="*70)
    print("\nThis script compares:")
    print("  1. Approximate: ρ_approx(r) = φ_HOMO(r) × φ_LUMO(r)")
    print("  2. Exact (TDDFT): ρ_exact(r) from transition density matrix")
    print("  3. Analytical: ρ_analytical(r) from HOMO⊗LUMO matrix")
    
    # Check if files exist
    print(f"\n{'='*70}")
    print("CHECKING INPUT FILES")
    print(f"{'='*70}")
    
    required_files = [HOMO_CUBE, LUMO_CUBE, TDDFT_TRANSITION_CUBE]
    optional_files = [ANALYTICAL_CUBE]
    
    for f in required_files:
        if os.path.exists(f):
            print(f"  ✓ Found: {f}")
        else:
            print(f"  ✗ Missing: {f}")
            print(f"\nError: Required file not found!")
            print(f"Please run tdm_calc_accurate.py first with:")
            print(f"  GENERATE_HOMO_LUMO = True")
            print(f"  GENERATE_TRANSITION_DENSITY = True")
            print(f"  STATES_TO_OUTPUT = [0]")
            sys.exit(1)
    
    has_analytical = os.path.exists(ANALYTICAL_CUBE)
    if has_analytical:
        print(f"  ✓ Found: {ANALYTICAL_CUBE}")
    else:
        print(f"  ⚠ Not found: {ANALYTICAL_CUBE}")
        print(f"    (This is optional, will skip analytical comparison)")
    
    # Read cube files
    print(f"\n{'='*70}")
    print("READING CUBE FILES")
    print(f"{'='*70}")
    
    print("  Reading HOMO.cube...")
    atoms_homo, origin_homo, grid_homo, data_homo = read_cube_file(HOMO_CUBE)
    nx_homo, ny_homo, nz_homo = data_homo.shape
    print(f"    Grid: {nx_homo} × {ny_homo} × {nz_homo} = {nx_homo*ny_homo*nz_homo:,} points")
    
    print("  Reading LUMO.cube...")
    atoms_lumo, origin_lumo, grid_lumo, data_lumo = read_cube_file(LUMO_CUBE)
    nx_lumo, ny_lumo, nz_lumo = data_lumo.shape
    print(f"    Grid: {nx_lumo} × {ny_lumo} × {nz_lumo} = {nx_lumo*ny_lumo*nz_lumo:,} points")
    
    print("  Reading TDDFT transition_density_state1.cube...")
    atoms_tddft, origin_tddft, grid_tddft, data_tddft = read_cube_file(TDDFT_TRANSITION_CUBE)
    nx_tddft, ny_tddft, nz_tddft = data_tddft.shape
    print(f"    Grid: {nx_tddft} × {ny_tddft} × {nz_tddft} = {nx_tddft*ny_tddft*nz_tddft:,} points")
    
    if has_analytical:
        print("  Reading analytical transition_HOMO_LUMO_analytical.cube...")
        atoms_anal, origin_anal, grid_anal, data_anal = read_cube_file(ANALYTICAL_CUBE)
        nx_anal, ny_anal, nz_anal = data_anal.shape
        print(f"    Grid: {nx_anal} × {ny_anal} × {nz_anal} = {nx_anal*ny_anal*nz_anal:,} points")
    
    # Verify grids match
    print(f"\n{'='*70}")
    print("VERIFYING GRID COMPATIBILITY")
    print(f"{'='*70}")
    
    if (nx_homo, ny_homo, nz_homo) != (nx_lumo, ny_lumo, nz_lumo):
        print("  ✗ Error: HOMO and LUMO grids don't match!")
        sys.exit(1)
    
    if (nx_homo, ny_homo, nz_homo) != (nx_tddft, ny_tddft, nz_tddft):
        print("  ✗ Error: HOMO/LUMO and TDDFT grids don't match!")
        sys.exit(1)
    
    if has_analytical and (nx_homo, ny_homo, nz_homo) != (nx_anal, ny_anal, nz_anal):
        print("  ✗ Error: Analytical grid doesn't match!")
        sys.exit(1)
    
    print(f"  ✓ All grids match: {nx_homo} × {ny_homo} × {nz_homo}")
    
    # Perform multiplication
    print(f"\n{'='*70}")
    print("CALCULATING APPROXIMATE TRANSITION DENSITY")
    print(f"{'='*70}")
    print("  Computing: ρ_approx(r) = φ_HOMO(r) × φ_LUMO(r)")
    
    # Simple point-wise multiplication
    data_multiplied = data_homo * data_lumo
    
    # Note: For transition density, we should actually compute:
    # ρ_trans(r) = φ_HOMO(r) × φ_LUMO(r) + φ_LUMO(r) × φ_HOMO(r) = 2 × φ_HOMO(r) × φ_LUMO(r)
    # But since multiplication is commutative, we just multiply by 2
    data_multiplied = 2.0 * data_multiplied
    
    print(f"  ✓ Multiplication complete")
    print(f"    Result shape: {data_multiplied.shape}")
    print(f"    Min value: {np.min(data_multiplied):12.6e}")
    print(f"    Max value: {np.max(data_multiplied):12.6e}")
    
    # Write output cube file
    print(f"\n{'='*70}")
    print("WRITING OUTPUT CUBE FILE")
    print(f"{'='*70}")
    
    comment1 = "Approximate transition density: 2 × φ_HOMO(r) × φ_LUMO(r)"
    comment2 = "Generated by compare_homo_lumo_multiplication.py"
    
    write_cube_file(OUTPUT_CUBE, atoms_homo, origin_homo, grid_homo, 
                    data_multiplied, comment1, comment2)
    print(f"  ✓ Written: {OUTPUT_CUBE}")
    
    # Compare with TDDFT result
    print(f"\n{'='*70}")
    print("COMPARISON 1: MULTIPLIED vs TDDFT")
    print(f"{'='*70}")
    
    stats1 = calculate_statistics(data_multiplied, data_tddft, 
                                   "Multiplied (approx)", "TDDFT (exact)")
    
    # Compare with analytical result (if available)
    if has_analytical:
        print(f"\n{'='*70}")
        print("COMPARISON 2: MULTIPLIED vs ANALYTICAL")
        print(f"{'='*70}")
        
        stats2 = calculate_statistics(data_multiplied, data_anal,
                                       "Multiplied (approx)", "Analytical (exact)")
        
        print(f"\n{'='*70}")
        print("COMPARISON 3: TDDFT vs ANALYTICAL")
        print(f"{'='*70}")
        print("(These should be nearly identical)")
        
        stats3 = calculate_statistics(data_tddft, data_anal,
                                       "TDDFT (exact)", "Analytical (exact)")
    
    # Generate difference cube files
    if GENERATE_DIFFERENCE_CUBES:
        print(f"\n{'='*70}")
        print("GENERATING DIFFERENCE CUBE FILES")
        print(f"{'='*70}")
        
        # Difference: Multiplied - TDDFT
        diff1 = data_multiplied - data_tddft
        diff1_file = os.path.join(CUBE_DIR, 'difference_multiplied_minus_tddft.cube')
        write_cube_file(diff1_file, atoms_homo, origin_homo, grid_homo, diff1,
                       "Difference: Multiplied - TDDFT",
                       "Shows error in approximation ρ ≈ 2×φ_HOMO×φ_LUMO")
        print(f"  ✓ Written: {diff1_file}")
        
        if has_analytical:
            # Difference: Multiplied - Analytical
            diff2 = data_multiplied - data_anal
            diff2_file = os.path.join(CUBE_DIR, 'difference_multiplied_minus_analytical.cube')
            write_cube_file(diff2_file, atoms_homo, origin_homo, grid_homo, diff2,
                           "Difference: Multiplied - Analytical",
                           "Shows error in approximation vs analytical HOMO⊗LUMO")
            print(f"  ✓ Written: {diff2_file}")
            
            # Difference: TDDFT - Analytical (should be very small)
            diff3 = data_tddft - data_anal
            diff3_file = os.path.join(CUBE_DIR, 'difference_tddft_minus_analytical.cube')
            write_cube_file(diff3_file, atoms_homo, origin_homo, grid_homo, diff3,
                           "Difference: TDDFT - Analytical",
                           "Should be nearly zero (both are exact methods)")
            print(f"  ✓ Written: {diff3_file}")
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    print("\nGenerated files:")
    print(f"  1. {OUTPUT_CUBE}")
    print(f"     Approximate transition density from HOMO×LUMO multiplication")
    
    if GENERATE_DIFFERENCE_CUBES:
        print(f"\n  2. {os.path.join(CUBE_DIR, 'difference_multiplied_minus_tddft.cube')}")
        print(f"     Difference showing approximation error")
        
        if has_analytical:
            print(f"\n  3. {os.path.join(CUBE_DIR, 'difference_multiplied_minus_analytical.cube')}")
            print(f"     Difference vs analytical method")
            print(f"\n  4. {os.path.join(CUBE_DIR, 'difference_tddft_minus_analytical.cube')}")
            print(f"     Difference between two exact methods (should be ~0)")
    
    print("\nKey findings:")
    print(f"  Similarity (Multiplied vs TDDFT): {stats1['similarity']:.6f}")
    
    if stats1['similarity'] > 0.95:
        print(f"  ✓ Approximation is EXCELLENT (similarity > 0.95)")
    elif stats1['similarity'] > 0.85:
        print(f"  ✓ Approximation is GOOD (similarity > 0.85)")
    elif stats1['similarity'] > 0.70:
        print(f"  ⚠ Approximation is FAIR (similarity > 0.70)")
    else:
        print(f"  ⚠ Approximation is POOR (similarity < 0.70)")
    
    print("\nInterpretation:")
    print("  - High similarity (>0.9): The approximation φ_HOMO×φ_LUMO works well")
    print("  - Lower similarity (<0.9): Basis function overlap effects are significant")
    print("  - The exact methods (TDDFT and Analytical) should be nearly identical")
    
    print("\nVisualization:")
    print("  Load all cube files in VMD/Jmol to see the differences visually:")
    print(f"    vmd {HOMO_CUBE} {LUMO_CUBE} {OUTPUT_CUBE} {TDDFT_TRANSITION_CUBE}")
    
    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
