# Guide: Charged and Open-Shell Calculations

## Overview

The script now supports calculations for:
- **Neutral molecules** (charge = 0)
- **Cations** (charge = +1, +2, ...)
- **Anions** (charge = -1, -2, ...)
- **Closed-shell systems** (spin = 1, singlet)
- **Open-shell systems** (spin = 2, doublet; spin = 3, triplet, etc.)

---

## Configuration

### Basic Settings

```python
# --- Charge and Spin Settings ---
CHARGE = 0  # Molecular charge: 0 (neutral), +1 (cation), -1 (anion)
SPIN = None  # Spin multiplicity (2S+1): None = auto-calculate
```

### Examples

#### Neutral Molecule (Default)
```python
CHARGE = 0
SPIN = None  # Auto-calculates: 1 (singlet) for even electrons
```

#### Cation (+1)
```python
CHARGE = +1
SPIN = None  # Auto-calculates: 2 (doublet) for odd electrons
```

#### Anion (-1)
```python
CHARGE = -1
SPIN = None  # Auto-calculates: 2 (doublet) for odd electrons
```

#### Triplet State
```python
CHARGE = 0
SPIN = 3  # Triplet (2 unpaired electrons)
```

---

## How It Works

### 1. Automatic Spin Calculation

If `SPIN = None`, the script automatically calculates spin multiplicity:

```
Number of electrons = N_neutral - charge

If even electrons → spin = 1 (singlet, closed-shell)
If odd electrons → spin = 2 (doublet, open-shell)
```

**Example for PTCDA:**
- Neutral: 166 electrons → spin = 1 (singlet)
- Cation (+1): 165 electrons → spin = 2 (doublet)
- Anion (-1): 167 electrons → spin = 2 (doublet)

### 2. DFT Method Selection

The script automatically selects the appropriate DFT method:

| Spin | System Type | DFT Method |
|------|-------------|------------|
| 1 | Closed-shell (singlet) | RKS (Restricted Kohn-Sham) |
| >1 | Open-shell (doublet, triplet, ...) | UKS (Unrestricted Kohn-Sham) |

### 3. TDDFT Calculation

TDDFT automatically uses the correct method based on the ground state:
- **RKS → TDDFT (closed-shell)**
- **UKS → TDDFT (open-shell)**

---

## Output Information

### Molecule Setup

```
======================================================================
MOLECULE SETUP
======================================================================
Loading molecule from: PTCDA.xyz
Basis set: 6-31g
Charge: 0
Number of atoms: 38
Number of electrons: 166
Number of basis functions: 302
Molecular charge: 0
Spin multiplicity (2S+1): 1
Number of unpaired electrons (2S): 0
System type: Closed-shell (singlet)
DFT method: RKS
```

**Key information:**
- **Molecular charge:** Confirms the charge setting
- **Spin multiplicity (2S+1):** 1 = singlet, 2 = doublet, 3 = triplet
- **Number of unpaired electrons (2S):** 0 = closed-shell, 1 = one unpaired, etc.
- **System type:** Closed-shell or open-shell
- **DFT method:** RKS (closed) or UKS (open)

### Verification

**Before calculation starts, verify:**
1. ✓ Charge is correct
2. ✓ Spin multiplicity is correct
3. ✓ Number of electrons is correct
4. ✓ DFT method is appropriate

**If incorrect, stop and adjust `CHARGE` or `SPIN` settings!**

---

## Examples

### Example 1: Neutral PTCDA (Singlet)

```python
CHARGE = 0
SPIN = None  # Auto: spin = 1
```

**Output:**
```
Number of electrons: 166
Spin multiplicity: 1
System type: Closed-shell (singlet)
DFT method: RKS
```

**Use case:** Ground state absorption spectrum

### Example 2: PTCDA Cation (Doublet)

```python
CHARGE = +1
SPIN = None  # Auto: spin = 2
```

**Output:**
```
Number of electrons: 165
Spin multiplicity: 2
System type: Open-shell (doublet)
DFT method: UKS
```

**Use case:** Oxidized species, hole transport

### Example 3: PTCDA Anion (Doublet)

```python
CHARGE = -1
SPIN = None  # Auto: spin = 2
```

**Output:**
```
Number of electrons: 167
Spin multiplicity: 2
System type: Open-shell (doublet)
DFT method: UKS
```

**Use case:** Reduced species, electron transport

### Example 4: Neutral PTCDA Triplet

```python
CHARGE = 0
SPIN = 3  # Manual: triplet
```

**Output:**
```
Number of electrons: 166
Spin multiplicity: 3
System type: Open-shell (triplet)
DFT method: UKS
```

**Use case:** Triplet excited state, phosphorescence

---

## Technical Details

### Spin Multiplicity

**Definition:** Spin multiplicity = 2S + 1

where S = total spin angular momentum

| Multiplicity | Name | S | Unpaired e⁻ |
|--------------|------|---|-------------|
| 1 | Singlet | 0 | 0 |
| 2 | Doublet | 1/2 | 1 |
| 3 | Triplet | 1 | 2 |
| 4 | Quartet | 3/2 | 3 |
| 5 | Quintet | 2 | 4 |

### PySCF Convention

**Important:** PySCF uses `spin = 2S` (number of unpaired electrons), not `2S+1`!

The script handles this conversion automatically:
```python
mol.spin = SPIN - 1  # Convert multiplicity to PySCF format
```

### RKS vs UKS

**RKS (Restricted Kohn-Sham):**
- For closed-shell systems (all electrons paired)
- Single set of orbitals
- α and β electrons share same spatial orbitals

**UKS (Unrestricted Kohn-Sham):**
- For open-shell systems (unpaired electrons)
- Separate α and β orbitals
- α and β electrons have different spatial orbitals

---

## HOMO/LUMO for Open-Shell Systems

### Alpha vs Beta Orbitals

For **UKS calculations**, the script uses **alpha orbitals** for HOMO/LUMO:

```
Note: Using alpha orbitals for HOMO/LUMO (open-shell system)
HOMO index: 82 (alpha)
LUMO index: 83 (alpha)
```

**Why alpha?**
- Alpha orbitals typically contain the unpaired electron(s)
- HOMO→LUMO transitions usually involve alpha orbitals
- Consistent with most literature conventions

**Note:** For detailed analysis, you may want to examine both alpha and beta orbitals separately.

---

## Common Use Cases

### 1. Neutral Molecule Absorption

```python
CHARGE = 0
SPIN = None  # Singlet
```

**Calculates:** S₀ → S₁, S₂, ... transitions

### 2. Cation Absorption

```python
CHARGE = +1
SPIN = None  # Doublet
```

**Calculates:** D₀ → D₁, D₂, ... transitions (doublet states)

### 3. Anion Absorption

```python
CHARGE = -1
SPIN = None  # Doublet
```

**Calculates:** D₀ → D₁, D₂, ... transitions (doublet states)

### 4. Triplet-Triplet Absorption

```python
CHARGE = 0
SPIN = 3  # Triplet ground state
```

**Calculates:** T₁ → T₂, T₃, ... transitions

---

## Convergence Tips

### For Charged Systems

**If SCF doesn't converge:**

1. **Try level shifting:**
   ```python
   mf.level_shift = 0.5  # Add before mf.kernel()
   ```

2. **Use better initial guess:**
   ```python
   mf.init_guess = 'atom'  # or 'minao'
   ```

3. **Adjust DIIS:**
   ```python
   mf.diis_space = 12
   mf.diis_start_cycle = 1
   ```

### For Open-Shell Systems

**If UKS doesn't converge:**

1. **Use stability analysis:**
   ```python
   mf.stability()  # Check for instabilities
   ```

2. **Try different XC functional:**
   - B3LYP usually works well
   - Try PBE0 or ωB97X-D for difficult cases

3. **Use smaller basis first:**
   - Start with 6-31g
   - Then use 6-31g* or larger

---

## Limitations

### 1. Spin Contamination

For **UKS calculations**, spin contamination can occur:
- Check `⟨S²⟩` value in output
- Should be close to S(S+1)
- Large deviation indicates contamination

### 2. TDDFT for Open-Shell

**TDDFT for open-shell systems:**
- More complex than closed-shell
- May have convergence issues
- Results should be validated carefully

### 3. Triplet and Higher

**For triplet and higher multiplicities:**
- Automatic spin calculation assumes doublet for odd electrons
- **Must manually set `SPIN`** for triplet, quartet, etc.

---

## Validation

### Check Your Results

1. **Verify electron count:**
   ```
   N_electrons = N_neutral - charge
   ```

2. **Verify spin:**
   ```
   Even electrons → singlet (spin=1)
   Odd electrons → doublet (spin=2)
   ```

3. **Check convergence:**
   ```
   ✓ SCF converged
   ✓ TDDFT converged
   ```

4. **Examine orbital energies:**
   - HOMO-LUMO gap should be reasonable
   - Compare with literature if available

---

## Summary

### Configuration

```python
CHARGE = 0  # 0, +1, -1, +2, -2, ...
SPIN = None  # Auto-calculate (recommended)
```

### Automatic Behavior

1. **Calculates spin** from electron count
2. **Selects RKS/UKS** based on spin
3. **Runs TDDFT** with appropriate method
4. **Handles HOMO/LUMO** for both closed/open-shell

### Key Points

- ✅ **Neutral:** CHARGE=0, auto spin
- ✅ **Cation:** CHARGE=+1, auto spin
- ✅ **Anion:** CHARGE=-1, auto spin
- ✅ **Triplet:** CHARGE=0, SPIN=3 (manual)
- ✅ **Verify** charge and spin before calculation!

---

## Example Output

### Neutral PTCDA

```
Molecular charge: 0
Spin multiplicity: 1
System type: Closed-shell (singlet)
DFT method: RKS
Ground state energy: -1522.456789 a.u.
✓ SCF converged
TDDFT method: TDDFT (RKS-based)
✓ TDDFT converged
```

### PTCDA Cation

```
Molecular charge: 1
Spin multiplicity: 2
System type: Open-shell (doublet)
DFT method: UKS
Ground state energy: -1522.123456 a.u.
✓ SCF converged
TDDFT method: TDDFT (UKS-based)
Note: Using alpha orbitals for HOMO/LUMO (open-shell system)
✓ TDDFT converged
```

The script now fully supports charged and open-shell calculations! 🎉
