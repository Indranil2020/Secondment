# Frequently Asked Questions

## Q1: Does TDDFT need separate basis set and XC functional settings?

**Answer:** No, TDDFT automatically inherits both the basis set and XC functional from the ground state DFT calculation.

**How it works:**
```python
# Ground state DFT
mf = dft.RKS(mol)
mf.xc = 'b3lyp'      # Set XC functional
# Basis set comes from mol object

# TDDFT calculation
td = tddft.TDDFT(mf)  # Inherits everything from mf
```

The TDDFT object `td` automatically uses:
- The same basis set that was used to create the `mol` object
- The same XC functional that was set in `mf.xc`
- The same molecular orbitals from the ground state

**Now configurable:** The script now has `XC_FUNCTIONAL` setting at the top, so you can easily change it without modifying the code.

---

## Q2: Why is B3LYP hardcoded? Can I use other functionals?

**Answer:** It was hardcoded, but now it's configurable!

**New configuration option:**
```python
XC_FUNCTIONAL = 'b3lyp'  # Change this to any functional
```

**Recommended functionals for excited states:**

| Functional | Type | Best For | Notes |
|------------|------|----------|-------|
| `'b3lyp'` | Hybrid | General purpose | Good balance, widely used |
| `'pbe0'` | Hybrid | Excited states | Often better than B3LYP for TDDFT |
| `'cam-b3lyp'` | Range-separated | Charge transfer | Better for long-range excitations |
| `'wb97x-d'` | Range-separated | CT + dispersion | Includes dispersion corrections |
| `'pbe'` | GGA | Fast screening | Faster but less accurate |
| `'blyp'` | GGA | Fast calculations | Pure functional, no exact exchange |

**For PTCDA:**
- Start with `'b3lyp'` (standard)
- Try `'cam-b3lyp'` if you have charge-transfer states
- Try `'pbe0'` for comparison

**Example configurations:**

```python
# Standard (good for most cases)
XC_FUNCTIONAL = 'b3lyp'

# Better for charge transfer excitations
XC_FUNCTIONAL = 'cam-b3lyp'

# Alternative hybrid functional
XC_FUNCTIONAL = 'pbe0'

# With dispersion corrections
XC_FUNCTIONAL = 'wb97x-d'
```

---

## Q3: Are STATES_TO_OUTPUT and NTO_STATES the same?

**Answer:** No, they are intentionally different and serve different purposes!

### STATES_TO_OUTPUT
**Purpose:** Which states to generate **cube files** for

**Generates:**
- `transition_density_state*.cube`
- `excited_state_density_state*.cube`
- `density_difference_state*.cube`

**File size:** ~150-500 MB per state (large!)

**Example:**
```python
STATES_TO_OUTPUT = [0, 1, 2]  # Generate cube files for states 1, 2, 3
```

### NTO_STATES
**Purpose:** Which states to perform **NTO analysis** for

**Generates:**
- `nto_state_*.molden`

**File size:** ~5 MB per state (small!)

**Example:**
```python
NTO_STATES = [0, 1, 2, 3, 4]  # NTO analysis for states 1-5
```

---

## Why Have Two Separate Settings?

**Reason 1: File Size**
- Cube files are huge (~150-500 MB each)
- NTO molden files are tiny (~5 MB each)
- You might want NTO for many states but cube files for only a few

**Reason 2: Different Uses**
- **Cube files:** For visualizing electron density (VMD, Jmol)
- **NTO files:** For analyzing orbital contributions (Jmol, Avogadro)

**Reason 3: Computational Cost**
- Generating cube files is slow (grid evaluation)
- NTO analysis is fast (just orbital transformation)

---

## Common Usage Patterns

### Pattern 1: Minimal Output
```python
STATES_TO_OUTPUT = [0]        # Only first state cube files
NTO_STATES = [0, 1, 2]        # NTO for first 3 states
```
**Result:** 3 cube files (~450 MB) + 3 molden files (~15 MB)

### Pattern 2: Detailed Analysis
```python
STATES_TO_OUTPUT = [0, 1, 2]  # Cube files for first 3 states
NTO_STATES = range(10)        # NTO for all 10 states
```
**Result:** 9 cube files (~1.4 GB) + 10 molden files (~50 MB)

### Pattern 3: Focus on Specific States
```python
STATES_TO_OUTPUT = [0, 4]     # Only states 1 and 5
NTO_STATES = [0, 1, 2, 3, 4]  # NTO for states 1-5
```
**Result:** 6 cube files (~900 MB) + 5 molden files (~25 MB)

### Pattern 4: Only NTO Analysis
```python
STATES_TO_OUTPUT = []         # No cube files
NTO_STATES = range(10)        # NTO for all states
GENERATE_TRANSITION_DENSITY = False
GENERATE_EXCITED_DENSITY = False
GENERATE_DENSITY_DIFFERENCE = False
```
**Result:** 0 cube files + 10 molden files (~50 MB)

---

## Quick Reference Table

| Setting | Purpose | Output | File Size | Typical Use |
|---------|---------|--------|-----------|-------------|
| `STATES_TO_OUTPUT` | Cube file generation | `.cube` files | 150-500 MB/state | Density visualization |
| `NTO_STATES` | NTO analysis | `.molden` files | ~5 MB/state | Orbital analysis |

---

## Recommendations

### For PTCDA:

**Quick test:**
```python
XC_FUNCTIONAL = 'b3lyp'
STATES_TO_OUTPUT = [0]
NTO_STATES = [0, 1, 2]
```

**Standard calculation:**
```python
XC_FUNCTIONAL = 'b3lyp'
STATES_TO_OUTPUT = [0, 1, 2]
NTO_STATES = [0, 1, 2, 3, 4]
```

**High quality:**
```python
XC_FUNCTIONAL = 'cam-b3lyp'  # Better for charge transfer
STATES_TO_OUTPUT = [0, 1, 2, 3, 4]
NTO_STATES = range(10)
```

---

## Additional Tips

### Choosing XC Functional

1. **Start with B3LYP** - Standard choice, well-tested
2. **Use CAM-B3LYP** if:
   - You have charge-transfer excitations
   - Excitations involve spatially separated regions
   - B3LYP gives unrealistic low energies

3. **Use PBE0** if:
   - You want an alternative hybrid functional
   - Literature uses PBE0 for comparison

4. **Use ωB97X-D** if:
   - You need dispersion corrections
   - System has π-π stacking

### Choosing States

1. **STATES_TO_OUTPUT:**
   - Start with `[0]` (first state only)
   - Add more if you need to visualize multiple states
   - Remember: each state = ~150-500 MB

2. **NTO_STATES:**
   - Include all states you want to understand
   - NTO files are small, so be generous
   - Useful for identifying orbital contributions

### Disk Space Planning

| Configuration | Cube Files | NTO Files | Total |
|---------------|------------|-----------|-------|
| Minimal (1 state) | ~450 MB | ~15 MB | ~465 MB |
| Standard (3 states) | ~1.4 GB | ~25 MB | ~1.4 GB |
| Full (5 states + 10 NTO) | ~2.3 GB | ~50 MB | ~2.4 GB |

---

## Summary

1. ✅ **TDDFT inherits basis set and XC functional** from ground state DFT
2. ✅ **XC functional is now configurable** - change `XC_FUNCTIONAL` at the top
3. ✅ **STATES_TO_OUTPUT and NTO_STATES are different:**
   - `STATES_TO_OUTPUT`: Cube files (large, for density visualization)
   - `NTO_STATES`: Molden files (small, for orbital analysis)
4. ✅ **You can have different lists** for each - this is intentional and useful!

All settings are clearly documented in the configuration section with comments explaining the differences.
