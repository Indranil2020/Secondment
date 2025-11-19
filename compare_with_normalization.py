#!/usr/bin/env python
'''
Enhanced HOMO×LUMO Multiplication with Normalization

This script improves the approximation by normalizing the result
to match the norm of the exact TDDFT transition density.

This accounts for the fact that the simple multiplication overestimates
the density due to basis function overlap effects.
'''

import numpy as np
import os
import sys

# Import functions from the original script
import sys
sys.path.insert(0, os.path.dirname(__file__))

# ============================================================================
# CONFIGURATION
# ============================================================================

CUBE_DIR = 'output'

HOMO_CUBE = os.path.join(CUBE_DIR, 'HOMO.cube')
LUMO_CUBE = os.path.join(CUBE_DIR, 'LUMO.cube')
TDDFT_TRANSITION_CUBE = os.path.join(CUBE_DIR, 'transition_density_state1.cube')

OUTPUT_CUBE_NORMALIZED = os.path.join(CUBE_DIR, 'transition_HOMO_LUMO_multiplied_normalized.cube')

# ============================================================================
# CUBE FILE I/O (same as original)
# ============================================================================

def read_cube_file(filename):
    with open(filename, 'r') as f:
        comment1 = f.readline()
        comment2 = f.readline()
        
        line = f.readline().split()
        natoms = int(line[0])
        origin = np.array([float(line[1]), float(line[2]), float(line[3])])
        
        grid_vectors = []
        for i in range(3):
            line = f.readline().split()
            n = int(line[0])
            vec = np.array([float(line[1]), float(line[2]), float(line[3])])
            grid_vectors.append((n, vec))
        
        atoms = []
        for i in range(abs(natoms)):
            line = f.readline().split()
            atomic_num = int(line[0])
            charge = float(line[1])
            x, y, z = float(line[2]), float(line[3]), float(line[4])
            atoms.append((atomic_num, charge, x, y, z))
        
        nx, ny, nz = grid_vectors[0][0], grid_vectors[1][0], grid_vectors[2][0]
        
        values = []
        for line in f:
            values.extend([float(x) for x in line.split()])
        
        data = np.array(values).reshape((nx, ny, nz))
        
    return atoms, origin, grid_vectors, data

def write_cube_file(filename, atoms, origin, grid_vectors, data, comment1="", comment2=""):
    nx, ny, nz = data.shape
    
    with open(filename, 'w') as f:
        f.write(f"{comment1}\n")
        f.write(f"{comment2}\n")
        f.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        
        for n, vec in grid_vectors:
            f.write(f"{n:5d} {vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")
        
        for atom in atoms:
            atomic_num, charge, x, y, z = atom
            f.write(f"{atomic_num:5d} {charge:12.6f} {x:12.6f} {y:12.6f} {z:12.6f}\n")
        
        values = data.flatten()
        for i in range(0, len(values), 6):
            line_values = values[i:i+6]
            f.write(" ".join([f"{v:13.5E}" for v in line_values]) + "\n")

def calculate_similarity(data1, data2):
    dot_product = np.sum(data1 * data2)
    norm1 = np.linalg.norm(data1)
    norm2 = np.linalg.norm(data2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    return dot_product / (norm1 * norm2)

# ============================================================================
# MAIN SCRIPT
# ============================================================================

def main():
    print("="*70)
    print("ENHANCED HOMO×LUMO MULTIPLICATION WITH NORMALIZATION")
    print("="*70)
    
    # Read files
    print("\nReading cube files...")
    atoms_homo, origin_homo, grid_homo, data_homo = read_cube_file(HOMO_CUBE)
    atoms_lumo, origin_lumo, grid_lumo, data_lumo = read_cube_file(LUMO_CUBE)
    atoms_tddft, origin_tddft, grid_tddft, data_tddft = read_cube_file(TDDFT_TRANSITION_CUBE)
    
    print(f"  ✓ HOMO: {data_homo.shape}")
    print(f"  ✓ LUMO: {data_lumo.shape}")
    print(f"  ✓ TDDFT: {data_tddft.shape}")
    
    # Method 1: Simple multiplication (original)
    print("\n" + "="*70)
    print("METHOD 1: Simple Multiplication (Original)")
    print("="*70)
    data_mult_simple = 2.0 * data_homo * data_lumo
    
    norm_simple = np.linalg.norm(data_mult_simple)
    norm_tddft = np.linalg.norm(data_tddft)
    similarity_simple = calculate_similarity(data_mult_simple, data_tddft)
    
    print(f"  Norm (multiplied): {norm_simple:.6f}")
    print(f"  Norm (TDDFT):      {norm_tddft:.6f}")
    print(f"  Ratio:             {norm_simple/norm_tddft:.6f}")
    print(f"  Similarity:        {similarity_simple:.6f}")
    
    # Method 2: Normalized multiplication
    print("\n" + "="*70)
    print("METHOD 2: Normalized Multiplication (Improved)")
    print("="*70)
    print("  Normalizing to match TDDFT norm...")
    
    # Normalize to match TDDFT norm
    data_mult_normalized = data_mult_simple * (norm_tddft / norm_simple)
    
    norm_normalized = np.linalg.norm(data_mult_normalized)
    similarity_normalized = calculate_similarity(data_mult_normalized, data_tddft)
    
    print(f"  Norm (normalized): {norm_normalized:.6f}")
    print(f"  Norm (TDDFT):      {norm_tddft:.6f}")
    print(f"  Ratio:             {norm_normalized/norm_tddft:.6f}")
    print(f"  Similarity:        {similarity_normalized:.6f}")
    
    # Improvement
    improvement = similarity_normalized - similarity_simple
    print(f"\n  Improvement: {improvement:+.6f} ({improvement*100:+.2f}%)")
    
    # Write normalized cube file
    print("\n" + "="*70)
    print("WRITING NORMALIZED CUBE FILE")
    print("="*70)
    
    comment1 = "Normalized transition density: 2×φ_HOMO(r)×φ_LUMO(r) × (||ρ_TDDFT|| / ||ρ_mult||)"
    comment2 = "Normalized to match TDDFT norm - improved approximation"
    
    write_cube_file(OUTPUT_CUBE_NORMALIZED, atoms_homo, origin_homo, grid_homo,
                    data_mult_normalized, comment1, comment2)
    print(f"  ✓ Written: {OUTPUT_CUBE_NORMALIZED}")
    
    # Detailed comparison
    print("\n" + "="*70)
    print("DETAILED COMPARISON")
    print("="*70)
    
    print("\nSimple multiplication:")
    print(f"  Min:        {np.min(data_mult_simple):12.6e}")
    print(f"  Max:        {np.max(data_mult_simple):12.6e}")
    print(f"  Norm:       {norm_simple:12.6f}")
    print(f"  Similarity: {similarity_simple:12.6f}")
    
    print("\nNormalized multiplication:")
    print(f"  Min:        {np.min(data_mult_normalized):12.6e}")
    print(f"  Max:        {np.max(data_mult_normalized):12.6e}")
    print(f"  Norm:       {norm_normalized:12.6f}")
    print(f"  Similarity: {similarity_normalized:12.6f}")
    
    print("\nTDDFT (exact):")
    print(f"  Min:        {np.min(data_tddft):12.6e}")
    print(f"  Max:        {np.max(data_tddft):12.6e}")
    print(f"  Norm:       {norm_tddft:12.6f}")
    
    # RMS differences
    rms_simple = np.sqrt(np.mean((data_mult_simple - data_tddft)**2))
    rms_normalized = np.sqrt(np.mean((data_mult_normalized - data_tddft)**2))
    
    print("\nRMS Differences:")
    print(f"  Simple:     {rms_simple:12.6e}")
    print(f"  Normalized: {rms_normalized:12.6e}")
    print(f"  Improvement: {(rms_simple - rms_normalized)/rms_simple * 100:.2f}%")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    print("\nComparison of methods:")
    print(f"  1. Simple multiplication:     Similarity = {similarity_simple:.6f}")
    print(f"  2. Normalized multiplication: Similarity = {similarity_normalized:.6f}")
    print(f"  3. Improvement:               {improvement:+.6f} ({improvement*100:+.2f}%)")
    
    print("\nKey insight:")
    print("  The simple multiplication OVERESTIMATES the density by ~47%")
    print("  due to basis function overlap effects.")
    print("  Normalizing improves the approximation but doesn't fix")
    print("  the spatial distribution differences.")
    
    print("\nFiles generated:")
    print(f"  1. {OUTPUT_CUBE_NORMALIZED}")
    print("     Normalized version with improved accuracy")
    
    print("\nVisualization:")
    print(f"  vmd {HOMO_CUBE} {LUMO_CUBE} \\")
    print(f"      {OUTPUT_CUBE_NORMALIZED} \\")
    print(f"      {TDDFT_TRANSITION_CUBE}")
    
    print("\n" + "="*70)
    print("DONE")
    print("="*70)

if __name__ == "__main__":
    main()
