# AO Basis vs Real Space: Clarifying the Confusion

## The Confusion

There are **two different representations** of the transition density:

1. **AO Basis Representation** (matrix): T_μν
2. **Real Space Representation** (cube file): ρ_trans(r)

They are related but **NOT the same thing**.

---

## Part 1: AO Basis Representation (Matrix)

### What is it?

The transition density matrix in the **atomic orbital (AO) basis**:

```
T_μν = C_μ^HOMO × C_ν^LUMO + C_μ^LUMO × C_ν^HOMO
```

### Properties:

- **Type:** Matrix (2D array)
- **Size:** N_AO × N_AO (e.g., 302 × 302 for PTCDA with 6-31g)
- **Elements:** Numbers representing density between basis functions μ and ν
- **Construction:** Outer product of MO coefficient vectors

### Example:

```python
# HOMO coefficients (vector)
C_HOMO = [0.1, 0.2, -0.3, 0.4, ...]  # Length = N_AO

# LUMO coefficients (vector)  
C_LUMO = [0.05, -0.15, 0.25, 0.1, ...]  # Length = N_AO

# Outer product creates matrix
T_μν = np.outer(C_HOMO, C_LUMO) + np.outer(C_LUMO, C_HOMO)
# Result: N_AO × N_AO matrix

# Example elements:
T[0,0] = C_HOMO[0] × C_LUMO[0] + C_LUMO[0] × C_HOMO[0] = 2 × 0.1 × 0.05 = 0.01
T[0,1] = C_HOMO[0] × C_LUMO[1] + C_LUMO[0] × C_HOMO[1] = 0.1×(-0.15) + 0.05×0.2 = -0.005
...
```

### What does T_μν mean?

**T_μν represents the transition density between basis functions μ and ν.**

- If μ and ν are on the same atom → local density change
- If μ and ν are on different atoms → non-local density correlation

**This is NOT a 3D spatial distribution yet!** It's a matrix in basis function space.

---

## Part 2: Real Space Representation (Cube File)

### What is it?

The transition density evaluated at each point in 3D space:

```
ρ_trans(r) = Σ_μν T_μν × χ_μ(r) × χ_ν(r)
```

where:
- **r** = (x, y, z) position in 3D space
- **χ_μ(r)** = atomic basis function μ evaluated at position r
- **T_μν** = transition density matrix element

### Properties:

- **Type:** 3D scalar field
- **Size:** nx × ny × nz grid points (e.g., 80 × 80 × 80 = 512,000 points)
- **Elements:** Density values at each spatial position
- **Construction:** Sum over all basis function pairs, weighted by T_μν

### Example:

```python
# At grid point r = (x=1.5, y=2.0, z=0.5) Angstrom

# Evaluate all basis functions at this point
χ_0(r) = 0.02  # Basis function 0 at position r
χ_1(r) = 0.15  # Basis function 1 at position r
χ_2(r) = -0.08 # Basis function 2 at position r
...

# Calculate density at this point
ρ_trans(r) = T[0,0]×χ_0(r)×χ_0(r) + T[0,1]×χ_0(r)×χ_1(r) + T[0,2]×χ_0(r)×χ_2(r) + ...
           = 0.01×0.02×0.02 + (-0.005)×0.02×0.15 + ...
           = (sum over all N_AO² terms)

# Repeat for all grid points (80×80×80 = 512,000 times)
```

### What does ρ_trans(r) mean?

**ρ_trans(r) is the transition density at position r in real space.**

- Positive values → electron accumulation (to LUMO)
- Negative values → electron depletion (from HOMO)
- This is what you visualize in VMD/Jmol

---

## Part 3: The Connection

### From Matrix to Cube File

```
AO Basis (Matrix)              Real Space (Cube File)
     T_μν                             ρ_trans(r)
  (N_AO × N_AO)                    (nx × ny × nz)
       ↓                                  ↑
       └──────── Basis function ─────────┘
                 evaluation
```

**The transformation:**

```python
# For each grid point r:
ρ_trans(r) = Σ_μ Σ_ν T_μν × χ_μ(r) × χ_ν(r)
```

This is a **double sum** over all basis functions.

### What PySCF Does

```python
cubegen.density(mol, filename, T_matrix, nx=nx, ny=ny, nz=nz)
```

**Internally:**
1. Creates 3D grid (nx × ny × nz points)
2. For each grid point r:
   - Evaluates all basis functions χ_μ(r)
   - Computes ρ(r) = Σ_μν T_μν χ_μ(r) χ_ν(r)
3. Writes to cube file

**This is computationally expensive!** For 302 basis functions and 512,000 grid points:
- 302² × 512,000 ≈ 46 billion operations

---

## Part 4: The Simplified Picture (Approximation)

### For Molecular Orbitals

When you generate a cube file for a **single orbital** (like HOMO):

```python
cubegen.orbital(mol, filename, C_HOMO, nx=nx, ny=ny, nz=nz)
```

This computes:
```
φ_HOMO(r) = Σ_μ C_μ^HOMO × χ_μ(r)
```

### For Transition Density (Simplified)

**Qualitatively**, you can think of:

```
ρ_trans(r) ≈ φ_HOMO(r) × φ_LUMO(r)
```

**But this is NOT exact!** The exact formula is:

```
ρ_trans(r) = Σ_μν T_μν × χ_μ(r) × χ_ν(r)
           = Σ_μν (C_μ^HOMO C_ν^LUMO + C_μ^LUMO C_ν^HOMO) × χ_μ(r) × χ_ν(r)
```

### Why the Approximation Works (Sort of)

If you expand the exact formula:

```
ρ_trans(r) = Σ_μν (C_μ^HOMO C_ν^LUMO + C_μ^LUMO C_ν^HOMO) × χ_μ(r) × χ_ν(r)
           = Σ_μ C_μ^HOMO χ_μ(r) × Σ_ν C_ν^LUMO χ_ν(r) + (symmetric term)
           = φ_HOMO(r) × φ_LUMO(r) + φ_LUMO(r) × φ_HOMO(r)
           = 2 × φ_HOMO(r) × φ_LUMO(r)
```

**Wait, this looks like it works!**

### The Catch

This derivation assumes the basis functions form a **complete orthonormal set**, which is only approximately true. In practice:

- Basis sets are **not orthogonal** (Gaussian basis functions overlap)
- Basis sets are **not complete** (finite basis set)
- The approximation is **qualitatively correct** but **quantitatively different**

---

## Part 5: Practical Implications

### What You Calculate in the Script

```python
# In AO basis (matrix)
homo_mo = mf.mo_coeff[:, homo_idx]  # Vector (N_AO,)
lumo_mo = mf.mo_coeff[:, lumo_idx]  # Vector (N_AO,)

# Create transition density matrix
t_homo_lumo = np.outer(homo_mo, lumo_mo) + np.outer(lumo_mo, homo_mo)
# Result: Matrix (N_AO, N_AO)

# Generate cube file (PySCF does the basis function evaluation)
cubegen.density(mol, filename, t_homo_lumo, nx=nx, ny=ny, nz=nz)
```

### What Gets Written to Cube File

The cube file contains:
- Grid dimensions (nx, ny, nz)
- Grid spacing (dx, dy, dz)
- Origin coordinates
- **ρ_trans(r) values at each grid point** (nx × ny × nz values)

### Verification

When you compare:
- `transition_density_state1.cube` (from TDDFT)
- `transition_HOMO_LUMO_analytical.cube` (from HOMO⊗LUMO)

You're comparing **real space densities** ρ(r), not matrices T_μν.

The similarity metric:
```python
similarity = Σ_r ρ_TDDFT(r) × ρ_analytical(r) / (||ρ_TDDFT|| × ||ρ_analytical||)
```

is computed on the **matrix representation** (before cube file generation) for efficiency:
```python
similarity = Σ_μν T_TDDFT[μ,ν] × T_analytical[μ,ν] / (||T_TDDFT|| × ||T_analytical||)
```

These give the same result because of the orthogonality of the basis function evaluation.

---

## Part 6: Summary Table

| Aspect | AO Basis (Matrix) | Real Space (Cube File) |
|--------|-------------------|------------------------|
| **Type** | Matrix | 3D scalar field |
| **Size** | N_AO × N_AO | nx × ny × nz |
| **Example** | 302 × 302 | 80 × 80 × 80 |
| **Elements** | Density between basis functions | Density at spatial points |
| **Construction** | Outer product of MO coefficients | Basis function evaluation |
| **Formula** | T_μν = C_μ^HOMO C_ν^LUMO + ... | ρ(r) = Σ_μν T_μν χ_μ(r) χ_ν(r) |
| **Used for** | Calculations, analysis | Visualization |
| **File format** | NumPy array (in memory) | Cube file (on disk) |

---

## Part 7: Common Misconceptions

### ❌ Misconception 1: "Cube file is just the matrix"

**Wrong!** The cube file is the matrix **evaluated in real space** using basis functions.

### ❌ Misconception 2: "ρ_trans(r) = φ_HOMO(r) × φ_LUMO(r) exactly"

**Wrong!** This is an approximation. The exact formula involves the full matrix contraction.

### ❌ Misconception 3: "I can multiply HOMO.cube and LUMO.cube to get transition density"

**Wrong!** You would get φ_HOMO(r) × φ_LUMO(r), which is close but not exact. The proper way is to construct T_μν first, then generate the cube file.

### ✅ Correct Understanding

1. Calculate T_μν in AO basis (matrix)
2. Use PySCF to evaluate ρ(r) = Σ_μν T_μν χ_μ(r) χ_ν(r)
3. Write ρ(r) to cube file
4. Visualize in VMD/Jmol

---

## Part 8: Why This Matters

### For Verification

When you verify S₁ is HOMO→LUMO:

**In AO basis:**
```python
T_TDDFT[μ,ν]  vs  T_analytical[μ,ν]
```

**In real space:**
```python
ρ_TDDFT(r)  vs  ρ_analytical(r)
```

Both should match if S₁ is pure HOMO→LUMO.

### For Understanding

- **Matrix T_μν**: Mathematical representation, used for calculations
- **Field ρ(r)**: Physical representation, used for visualization

You need both:
- Matrix for quantitative analysis (similarity, weights)
- Field for qualitative understanding (where electrons move)

---

## Part 9: Code Example

### Step 1: Create Matrix in AO Basis

```python
import numpy as np

# Get HOMO and LUMO MO coefficients (vectors)
homo_mo = mf.mo_coeff[:, homo_idx]  # Shape: (N_AO,)
lumo_mo = mf.mo_coeff[:, lumo_idx]  # Shape: (N_AO,)

# Create transition density matrix (outer product)
T_matrix = np.outer(homo_mo, lumo_mo) + np.outer(lumo_mo, homo_mo)
# Shape: (N_AO, N_AO)

print(f"Matrix shape: {T_matrix.shape}")  # e.g., (302, 302)
print(f"Matrix element T[0,0]: {T_matrix[0,0]}")
print(f"Matrix element T[0,1]: {T_matrix[0,1]}")
```

### Step 2: Generate Cube File (Real Space)

```python
from pyscf.tools import cubegen

# PySCF evaluates the matrix in real space
cubegen.density(mol, 'transition.cube', T_matrix, nx=80, ny=80, nz=80)

# This creates a cube file with 80×80×80 = 512,000 grid points
# Each point contains ρ(r) = Σ_μν T_matrix[μ,ν] × χ_μ(r) × χ_ν(r)
```

### Step 3: What's in the Cube File

```
# Cube file format (simplified)
HOMO-LUMO transition density
Generated by PySCF
38  -10.0  -10.0  -10.0    # N_atoms, origin
80   0.25   0.0   0.0      # nx, dx vector
80   0.0   0.25   0.0      # ny, dy vector  
80   0.0   0.0   0.25      # nz, dz vector
6  0.0  x1 y1 z1           # Atom 1 (C)
1  0.0  x2 y2 z2           # Atom 2 (H)
...
ρ(r1)  ρ(r2)  ρ(r3)  ...  # Density values at grid points
```

The density values are **ρ(r)**, not T_μν!

---

## Part 10: Final Clarification

### Your Diagram (Corrected)

```
AO Basis (coefficients)     →     Real Space (cube files)
     C_HOMO (vector)                φ_HOMO(r) (3D field)
     C_LUMO (vector)                φ_LUMO(r) (3D field)
           ↓                              ↓
    Outer product                  NOT simple multiplication!
    C_HOMO ⊗ C_LUMO               
         ↓                              
    T_μν (matrix)                  
         ↓                              
    Basis function evaluation
    Σ_μν T_μν χ_μ(r) χ_ν(r)
         ↓
    ρ_trans(r) (3D field)
```

### The Key Point

**You CANNOT go directly from φ_HOMO(r) and φ_LUMO(r) to ρ_trans(r) by simple multiplication.**

You must:
1. Go back to AO basis: construct T_μν from C_HOMO and C_LUMO
2. Evaluate in real space: ρ(r) = Σ_μν T_μν χ_μ(r) χ_ν(r)

**The script does this correctly!** When you call:
```python
cubegen.density(mol, filename, T_matrix, ...)
```

PySCF handles the basis function evaluation internally.

---

## Conclusion

- **AO basis (T_μν)**: Matrix representation, used for calculations
- **Real space (ρ(r))**: Field representation, used for visualization
- **Connection**: ρ(r) = Σ_μν T_μν χ_μ(r) χ_ν(r)
- **Approximation**: ρ(r) ≈ 2 φ_HOMO(r) φ_LUMO(r) (qualitatively correct)
- **Exact method**: Construct T_μν first, then evaluate with basis functions

The script implements the **exact method**, which is why the verification is reliable!
