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
import re
import os

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE SETTINGS
# ============================================================================

# --- Parallel Calculation Settings ---
ENABLE_PARALLEL = True  # Enable/disable parallel computation
NUM_THREADS = 0

# --- Molecule Selection ---
USE_XYZ = True
XYZ_FILE = 'PTCDA_clean.xyz'
BASIS_SET = '6-31g*'

# --- Charge and Spin Settings ---
CHARGE = -1
SPIN = None
# Note: Spin is auto-calculated from electron count if set to None
# For charged systems: cation (+1) typically has spin=2 (doublet), anion (-1) has spin=2 (doublet)
# Neutral even-electron systems typically have spin=1 (singlet)

ENABLE_CDFT = False
MONOMER_A_ATOMS = []
TARGET_CHARGE_A = 0.0
CDFT_VC_X0 = 0.0
CDFT_VC_X1 = 0.1
CDFT_CHARGE_TOL = 1e-4
CDFT_MAX_ITER = 25

# --- DFT/TDDFT Settings ---
# Note: TDDFT uses the same basis set and XC functional as ground state DFT
XC_FUNCTIONAL = 'wb97x'
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

# --- Calculation Stage Control ---
# These switches control which calculation stages to run
# Valid combinations:
#   OPTIMISE_GEOMETRY=True,  ENABLE_DFT=False, ENABLE_TDDFT=False  -> Optimization only
#   OPTIMISE_GEOMETRY=False, ENABLE_DFT=True,  ENABLE_TDDFT=False  -> DFT only (single point)
#   OPTIMISE_GEOMETRY=True,  ENABLE_DFT=True,  ENABLE_TDDFT=False  -> Optimization + DFT
#   OPTIMISE_GEOMETRY=True,  ENABLE_DFT=True,  ENABLE_TDDFT=True   -> Full workflow (Opt + DFT + TDDFT)
#   OPTIMISE_GEOMETRY=False, ENABLE_DFT=True,  ENABLE_TDDFT=True   -> DFT + TDDFT (no optimization)
# Note: TDDFT requires DFT (ENABLE_DFT is automatically set to True if ENABLE_TDDFT=True)
ENABLE_DFT = True               # Run ground state DFT calculation
ENABLE_TDDFT = True             # Run TDDFT excited state calculation (requires DFT)

# --- Emission Calculation Settings ---
# Emission = fluorescence from excited state minimum geometry
# Requires: 1) Optimize excited state geometry, 2) Calculate emission energy at that geometry
# ENABLE_EMISSION: Calculate emission energy (requires ENABLE_TDDFT=True)
# EMISSION_STATE: Which excited state to optimize for emission (0-indexed, typically 0 for S1)
# EMISSION_OPT_MAX_STEPS: Max steps for excited state geometry optimization
ENABLE_EMISSION = True
EMISSION_STATE = 0              # 0 = S1 (first excited state)
EMISSION_OPT_MAX_STEPS = 100    # Max steps for excited state optimization
EMISSION_OPT_CONV = 'tight'     # Convergence: 'tight', 'normal', 'loose'

# --- Geometry Optimization Settings ---
OPTIMISE_GEOMETRY = True
OPT_CYCLES = 1                  # Number of optimization cycles (output of cycle N feeds into cycle N+1)
OPT_MAX_STEPS = 100             # Max steps per optimization cycle
OPT_CONV_PARAMS = 'tight'       # Convergence: 'tight', 'normal', 'loose'

# --- Verbose/Debug Settings ---
VERBOSE_LEVEL = 4
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
GENERATE_DEFORMATION_DENSITY = True

# --- Output Directory ---
OUTPUT_DIR = 'PTCDA_clean_wb97x_6_31g__cpu_charge-1'

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
    
    WARNING: This assumes singlet for even electrons and doublet for odd.
    This may be incorrect for:
    - Diradicals (e.g., O2 is triplet despite even electrons)
    - Transition metal complexes with high-spin configurations
    - Open-shell singlets
    For such systems, explicitly set SPIN in the configuration.
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

margin_bohr = BOX_MARGIN / 0.529177

if USE_XYZ:
    print(f"Loading molecule from: {XYZ_FILE}")
    print(f"Basis set: {BASIS_SET}")
    print(f"Charge: {CHARGE}")
    mol, calculated_spin = create_molecule_from_xyz(XYZ_FILE, BASIS_SET, CHARGE, SPIN)
    actual_spin = calculated_spin
    # Warn if spin was auto-detected for even-electron system
    if SPIN is None and mol.nelectron % 2 == 0:
        print("  ⚠️  Note: Spin auto-detected as singlet for even-electron system.")
        print("     If this is a diradical/triplet (e.g., O2), set SPIN=3 explicitly.")
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

print(f"DFT method: {dft_method}")

# ==================================================================
# 2.  GEOMETRY  (optional) → GROUND STATE  →  PROPERTIES / TDDFT
# ==================================================================
print("\n" + "="*70)
print("GEOMETRY & GROUND STATE SETUP")
print("="*70)
# ---------- 2a.  parse functional / dispersion --------------------
def parse_xc_functional(xc_string):
    """
    Parse XC functional string to separate base functional and dispersion correction.
    
    Based on PySCF official source (pyscf/scf/dispersion.py):
    - Supported dispersion: d3bj, d3zero, d3bjm, d3zerom, d3op, d4
    - NOT supported: wb97x-d, wb97x-d3 (these are in PySCF's blacklist)
    - Special: wb97x-d3bj uses wb97x-v functional internally
    
    Examples:
    - 'b3lyp'           -> ('b3lyp', None)
    - 'b3lyp-d3bj'      -> ('b3lyp', 'd3bj')
    - 'wb97x-d3bj'      -> ('wb97x', 'd3bj')  # PySCF handles wb97x-v mapping
    - 'pbe0-d4'         -> ('pbe0', 'd4')
    - 'cam-b3lyp-d3bj'  -> ('cam-b3lyp', 'd3bj')
    
    Returns: (clean_xc, dispersion_type)
    """
    xc_lower = xc_string.lower().strip()
    
    # PySCF blacklist - these are NOT supported
    blacklist = {'wb97x-d', 'wb97x-d3', 'wb97x_d', 'wb97x_d3'}
    if xc_lower in blacklist:
        raise ValueError(
            f"'{xc_string}' is NOT supported by PySCF.\n"
            f"Use 'wb97x-d3bj' instead (recommended) or 'wb97x' without dispersion."
        )
    
    # Supported dispersion versions (from PySCF DISP_VERSIONS)
    # Order matters: longer/more specific patterns first
    dispersion_patterns = [
        (r'-d3bjm$', 'd3bjm'),         # d3bj modified
        (r'-d3bj$', 'd3bj'),           # d3bj (Becke-Johnson damping)
        (r'-d3zerom$', 'd3zerom'),     # d3zero modified  
        (r'-d3zero$', 'd3zero'),       # d3 zero-damping
        (r'-d3op$', 'd3op'),           # d3 optimized power
        (r'-d4$', 'd4'),               # d4
    ]
    
    for pattern, disp_type in dispersion_patterns:
        match = re.search(pattern, xc_lower)
        if match:
            clean_xc = xc_string[:match.start()]
            return clean_xc, disp_type
    
    # No dispersion suffix found
    return xc_string, None

clean_xc, disp_suffix = parse_xc_functional(XC_FUNCTIONAL)
print(f"Parsed XC: '{XC_FUNCTIONAL}' -> base='{clean_xc}', dispersion='{disp_suffix}'")

# ---------- 2b.  helper: build mf with current geometry ------------
def build_mf(mol, verbose_level=VERBOSE_LEVEL):
    """Build mf with identical settings for opt or final."""
    mf = (dft.RKS if actual_spin == 1 else dft.UKS)(mol, xc=clean_xc)
    if disp_suffix:
        mf.disp = disp_suffix
    # Set verbosity based on VERBOSE_LEVEL
    if verbose_level >= 3:
        mf.verbose = 5  # Debug: show everything
    elif verbose_level >= 2:
        mf.verbose = 4  # Detailed: show SCF iterations
    elif verbose_level >= 1:
        mf.verbose = 3  # Normal: show key info
    else:
        mf.verbose = 2  # Minimal
    return mf


def _cdft_dm_total(dm):
    if isinstance(dm, tuple):
        return dm[0] + dm[1]
    if hasattr(dm, 'ndim') and dm.ndim == 3 and dm.shape[0] == 2:
        return dm[0] + dm[1]
    return dm


def _cdft_build_weight_matrix(mol, atom_indices):
    atom_set = {int(x) for x in atom_indices}
    if len(atom_set) == 0:
        raise ValueError('ENABLE_CDFT=True but MONOMER_A_ATOMS is empty')

    ao_labels = mol.ao_labels(fmt=False)
    nao = mol.nao_nr()
    if len(ao_labels) != nao:
        raise ValueError(f'Unexpected ao_labels length {len(ao_labels)} != nao {nao}')

    ao_mask = np.fromiter((lbl[0] in atom_set for lbl in ao_labels), dtype=np.bool_, count=len(ao_labels))
    if not ao_mask.any():
        raise ValueError('MONOMER_A_ATOMS selects zero AOs (check 0-indexed atom numbering)')

    s1e = mol.intor('int1e_ovlp')
    w = 0.5 * (s1e * ao_mask[:, None] + s1e * ao_mask[None, :])

    neutral_elec = sum(mol.atom_charge(a) for a in atom_set)
    target_pop = float(neutral_elec) - float(TARGET_CHARGE_A)

    return w, target_pop, neutral_elec


def _cdft_population(w, dm):
    dm_tot = _cdft_dm_total(dm)
    dm_tot = np.asarray(dm_tot, dtype=np.float64)
    return np.einsum('ij,ji->', w, dm_tot).item()


def run_cdft_secant(mf, mol):
    print('Starting cDFT optimization loop...', flush=True)
    w, target_pop, neutral_elec = _cdft_build_weight_matrix(mol, MONOMER_A_ATOMS)
    print(f'Targeting exactly {target_pop:.6f} electrons on Monomer A (neutral = {neutral_elec}, target charge = {TARGET_CHARGE_A})', flush=True)

    orig_get_hcore = mf.get_hcore
    dm_prev = None

    def eval_diff(vc, dm0):
        print(f'  [cDFT] Testing Vc = {vc:.8f} ...', flush=True)

        def get_hcore_cdft(*args, **kwargs):
            hcore = np.asarray(orig_get_hcore(*args, **kwargs), dtype=np.float64)
            return hcore + (vc * w)

        mf.get_hcore = get_hcore_cdft
        mf.kernel(dm0=dm0)
        dm_new = mf.make_rdm1()
        pop_a = _cdft_population(w, dm_new)
        diff = pop_a - target_pop
        print(f'  [cDFT] Pop A = {pop_a:.6f} | Diff = {diff:.6e} | E = {mf.e_tot:.10f}', flush=True)
        return diff, dm_new

    vc0 = float(CDFT_VC_X0)
    vc1 = float(CDFT_VC_X1)
    f0, dm_prev = eval_diff(vc0, dm_prev)
    if abs(f0) < CDFT_CHARGE_TOL:
        mf.cdft_vc = vc0
        print(f'✓ cDFT converged (Vc = {vc0:.8f})', flush=True)
        return

    f1, dm_prev = eval_diff(vc1, dm_prev)
    if abs(f1) < CDFT_CHARGE_TOL:
        mf.cdft_vc = vc1
        print(f'✓ cDFT converged (Vc = {vc1:.8f})', flush=True)
        return

    for it in range(1, int(CDFT_MAX_ITER) + 1):
        denom = (f1 - f0)
        if abs(denom) < 1e-14:
            vc2 = vc1 + 0.1
        else:
            vc2 = vc1 - f1 * (vc1 - vc0) / denom

        vc0, f0 = vc1, f1
        vc1 = float(vc2)
        f1, dm_prev = eval_diff(vc1, dm_prev)

        if abs(f1) < CDFT_CHARGE_TOL:
            mf.cdft_vc = vc1
            print(f'✓ cDFT converged in {it} secant iterations', flush=True)
            print(f'✓ Optimal Vc multiplier: {vc1:.8f}', flush=True)
            return

    raise RuntimeError('cDFT optimization failed to reach the target charge constraint')

# ---------- 2c.  convergence parameters for optimization -----------
def get_opt_conv_params(conv_preset):
    """Get convergence parameters for geometry optimization."""
    presets = {
        'tight': {'convergence_energy': 1e-6, 'convergence_grms': 3e-4, 'convergence_gmax': 4.5e-4,
                  'convergence_drms': 1.2e-3, 'convergence_dmax': 1.8e-3},
        'normal': {'convergence_energy': 1e-5, 'convergence_grms': 1e-3, 'convergence_gmax': 1.5e-3,
                   'convergence_drms': 4e-3, 'convergence_dmax': 6e-3},
        'loose': {'convergence_energy': 1e-4, 'convergence_grms': 3e-3, 'convergence_gmax': 4.5e-3,
                  'convergence_drms': 1.2e-2, 'convergence_dmax': 1.8e-2},
    }
    if isinstance(conv_preset, dict):
        return conv_preset
    return presets.get(conv_preset.lower(), presets['normal'])

# ---------- 2d.  geometry optimisation? ---------------------------
if OPTIMISE_GEOMETRY:
    from pyscf.geomopt.geometric_solver import optimize as opt_kernel
    
    print("\n" + "-"*70)
    print("GEOMETRY OPTIMIZATION")
    print("-"*70)
    print(f"Optimization cycles: {OPT_CYCLES}")
    print(f"Max steps per cycle: {OPT_MAX_STEPS}")
    print(f"Convergence: {OPT_CONV_PARAMS}")
    if VERBOSE_LEVEL >= 2:
        conv_params = get_opt_conv_params(OPT_CONV_PARAMS)
        print(f"  Energy threshold:  {conv_params['convergence_energy']:.1e} a.u.")
        print(f"  Gradient RMS:      {conv_params['convergence_grms']:.1e} a.u./Bohr")
    print("-"*70)
    
    # Run multiple optimization cycles if requested
    for opt_cycle in range(1, OPT_CYCLES + 1):
        if OPT_CYCLES > 1:
            print(f"\n--- Optimization Cycle {opt_cycle}/{OPT_CYCLES} ---")
        
        print("Optimising geometry with dispersion-corrected forces...")
        mf = build_mf(mol)
        conv_params = get_opt_conv_params(OPT_CONV_PARAMS)
        mol = opt_kernel(mf, maxsteps=OPT_MAX_STEPS, **conv_params)
        
        # Save optimised structure (intermediate cycles get numbered)
        if OPT_CYCLES > 1:
            opt_xyz = os.path.join(OUTPUT_DIR, f'optimised_structure_cycle{opt_cycle}.xyz')
        else:
            opt_xyz = os.path.join(OUTPUT_DIR, 'optimised_structure.xyz')
        mol.tofile(opt_xyz, format='xyz')
        print(f"✓ Optimised geometry saved to: {opt_xyz}")
    
    # Always save final structure as optimised_structure.xyz for downstream use
    if OPT_CYCLES > 1:
        final_xyz = os.path.join(OUTPUT_DIR, 'optimised_structure.xyz')
        mol.tofile(final_xyz, format='xyz')
        print(f"✓ Final optimised geometry saved to: {final_xyz}")

# ---------- 2e.  Enforce calculation stage dependencies -----------------
# TDDFT requires DFT
if ENABLE_TDDFT and not ENABLE_DFT:
    print("Note: ENABLE_TDDFT=True requires DFT. Setting ENABLE_DFT=True automatically.")
    ENABLE_DFT = True

# If only optimization requested, skip DFT and TDDFT
if not ENABLE_DFT and not ENABLE_TDDFT:
    if OPTIMISE_GEOMETRY:
        print("\n" + "="*70)
        print("OPTIMIZATION ONLY MODE - Skipping DFT and TDDFT")
        print("="*70)
        print(f"Final optimised structure: {os.path.join(OUTPUT_DIR, 'optimised_structure.xyz')}")
        print("\n✓ Calculation completed (optimization only)")
        import sys
        sys.exit(0)
    else:
        print("ERROR: No calculation stages enabled. Set at least one of:")
        print("  - OPTIMISE_GEOMETRY=True")
        print("  - ENABLE_DFT=True")
        print("  - ENABLE_TDDFT=True")
        import sys
        sys.exit(1)

# ---------- 2f.  final SCF (only one kernel call) -----------------
if ENABLE_DFT:
    print("\n" + "-"*70)
    print("GROUND STATE DFT CALCULATION" + (" (WITH cDFT)" if ENABLE_CDFT else ""))
    print("-"*70)
    print(f"XC functional: {XC_FUNCTIONAL} (clean: {clean_xc})")
    print(f"Dispersion: {disp_suffix if disp_suffix else 'None'}")
    print(f"Basis set: {BASIS_SET}")
    print(f"Method: {dft_method}")
    print("-"*70)

    mf = build_mf(mol)

    if ENABLE_CDFT:
        run_cdft_secant(mf, mol)
    else:
        mf.kernel()

    if not mf.converged:
        print("⚠ WARNING: SCF did not converge!")
        print("  Try: 1) Different initial guess, 2) Level shifting, 3) DIIS settings")
    else:
        print("✓ SCF converged")

    # Print energy with precision based on verbose level
    if VERBOSE_LEVEL >= 2:
        print(f"✓ Ground-state energy: {mf.e_tot:.8f} a.u.")
    else:
        print(f"✓ Ground-state energy: {mf.e_tot:.6f} a.u.")

# ============================================================================
# 2. GROUND STATE DFT CALCULATION
# ============================================================================

# print("\n" + "="*70)
# print("GROUND STATE DFT CALCULATION")
# print("="*70)
# print(f"XC functional: {XC_FUNCTIONAL}")
# print(f"Basis set: {BASIS_SET}")
# print(f"Method: {dft_method}")

# # ---------- 1.  split “functional-d3bj” into (“functional”, “d3bj”) ---
# m = re.search(r'-((?:d3bj|d3|d3zero|d4|d3bps))$', XC_FUNCTIONAL.lower())
# clean_xc = XC_FUNCTIONAL[:m.start()] if m else XC_FUNCTIONAL   # remove suffix
# disp_suffix = m.group(1) if m else None                       # d3bj, d4, …

# # ---------- 2.  create RKS/UKS with the *clean* functional ------------
# if actual_spin == 1:
#     mf = dft.RKS(mol, xc=clean_xc)
# else:
#     mf = dft.UKS(mol, xc=clean_xc)

# # ---------- 3.  add dispersion if present -----------------------------
# if disp_suffix:
#     mf.disp = disp_suffix
#     print(f"✓ Adding DFT-{disp_suffix.upper()} dispersion correction")

# print('driver object :', type(mf))
# print('has disp attr :', hasattr(mf, 'disp'))
# # ---------- 4.  run SCF -------------------------------------------------
# mf.kernel()

# if not mf.converged:
#     print("WARNING: SCF did not converge!")
#     print("Try: 1) Different initial guess, 2) Level shifting, 3) DIIS settings")
# else:
#     print("✓ SCF converged")

# print(f"Ground state energy: {mf.e_tot:.6f} a.u.")

# ============================================================================
# 2A. GROUND STATE DENSITY AND ELECTROSTATIC POTENTIAL
# ============================================================================

if ENABLE_DFT and (GENERATE_GROUND_STATE_DENSITY or GENERATE_ELECTROSTATIC_POTENTIAL or GENERATE_DEFORMATION_DENSITY):
    print("\n" + "="*70)
    print("GROUND STATE DENSITY AND POTENTIAL")
    print("="*70)

    coords_ang = mol.atom_coords() * 0.529177
    mol_size_ang = coords_ang.max(axis=0) - coords_ang.min(axis=0)
    margin_bohr = BOX_MARGIN / 0.529177
    if USE_GRID_RESOLUTION:
        nx_g, ny_g, nz_g = GRID_RESOLUTION
    else:
        box_size_ang = mol_size_ang + 2 * BOX_MARGIN
        nx_g = int(np.ceil(box_size_ang[0] / GRID_SPACING))
        ny_g = int(np.ceil(box_size_ang[1] / GRID_SPACING))
        nz_g = int(np.ceil(box_size_ang[2] / GRID_SPACING))
    
    # Calculate ground state density matrix
    dm = mf.make_rdm1()

    # Overlap matrix for correct electron counting in AO basis
    s1e = mol.intor('int1e_ovlp')
    
    # Handle UKS (sum alpha and beta densities)
    # BOTH CPU and GPU PySCF return shape (2, nao, nao) for UKS!
    if isinstance(dm, tuple):
        # Old PySCF versions: tuple of (dm_alpha, dm_beta)
        dm_alpha, dm_beta = dm
        dm_alpha = np.asarray(dm_alpha, dtype=np.float64)
        dm_beta = np.asarray(dm_beta, dtype=np.float64)
        
        dm_total = dm_alpha + dm_beta
        
        n_alpha = np.einsum('ij,ji->', dm_alpha, s1e).item()
        n_beta = np.einsum('ij,ji->', dm_beta, s1e).item()
        
        print("System type: UKS (open-shell)")
        print(f"  Alpha electrons: {n_alpha:.2f}")
        print(f"  Beta electrons: {n_beta:.2f}")
    elif hasattr(dm, 'ndim') and dm.ndim == 3 and dm.shape[0] == 2:
        # Modern PySCF UKS: array with shape (2, nao, nao)
        dm = np.asarray(dm, dtype=np.float64)
        dm_alpha = dm[0]
        dm_beta = dm[1]
        
        dm_total = dm_alpha + dm_beta
        
        n_alpha = np.einsum('ij,ji->', dm_alpha, s1e).item()
        n_beta = np.einsum('ij,ji->', dm_beta, s1e).item()
        
        print("System type: UKS (open-shell)")
        print(f"  Alpha electrons: {n_alpha:.2f}")
        print(f"  Beta electrons: {n_beta:.2f}")
    else:
        # RKS: single 2D matrix
        dm_total = np.asarray(dm, dtype=np.float64)
        print("System type: RKS (closed-shell)")
    
    # Calculate total electrons
    total_electrons = np.einsum('ij,ji->', dm_total, s1e).item()
    print(f"Total electrons (Tr[DM*S]): {total_electrons:.6f}")
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
        cubegen.density(mol, density_file, dm_total, nx=nx_g, ny=ny_g, nz=nz_g, margin=margin_bohr)
        print(f"  ✓ Ground state density: {density_file}")
        print(f"    Use this to visualize total electron distribution")
    
    # Generate electrostatic potential cube file
    if GENERATE_ELECTROSTATIC_POTENTIAL:
        esp_file = os.path.join(OUTPUT_DIR, 'electrostatic_potential.cube')
        print(f"\nGenerating electrostatic potential (ESP)...")
        # ESP = Nuclear potential + Electronic potential
        # cubegen.mep calculates the molecular electrostatic potential
        cubegen.mep(mol, esp_file, dm_total, nx=nx_g, ny=ny_g, nz=nz_g, margin=margin_bohr)
        print(f"  ✓ Electrostatic potential: {esp_file}")
        print(f"    Use this to identify:")
        print(f"      - Nucleophilic sites (negative ESP, red)")
        print(f"      - Electrophilic sites (positive ESP, blue)")
        print(f"      - Reaction sites and molecular recognition")
    
    # Generate deformation density (SCF density - Promolecule density)
    if GENERATE_DEFORMATION_DENSITY:
        print(f"\nGenerating deformation density...")
        print(f"  Deformation density = SCF density - Promolecule density")
        print(f"  Shows charge redistribution due to chemical bonding")

        # Get promolecule density (superposition of atomic densities)
        # This is the non-interacting atomic density (no electron-electron interaction)
        if actual_spin == 1:
            # RKS: use scf.hf.init_guess_by_atom for closed-shell systems
            from pyscf import scf as pyscf_scf
            dm_promol = pyscf_scf.hf.init_guess_by_atom(mol)
            dm_promol = np.asarray(dm_promol, dtype=np.float64)

            promol_e = np.einsum('ij,ji->', dm_promol, s1e).item()
            print(f"  System: RKS (closed-shell)")
            print(f"  Promolecule electrons (Tr[DM*S]): {promol_e:.6f}")
        else:
            # UKS: use scf.uhf.init_guess_by_atom
            from pyscf import scf as pyscf_scf
            dm_promol = pyscf_scf.uhf.init_guess_by_atom(mol)

            # Handle both tuple and array formats
            if isinstance(dm_promol, tuple):
                dm_promol_alpha, dm_promol_beta = dm_promol
                dm_promol_alpha = np.asarray(dm_promol_alpha, dtype=np.float64)
                dm_promol_beta = np.asarray(dm_promol_beta, dtype=np.float64)
                dm_promol_total = dm_promol_alpha + dm_promol_beta
            elif dm_promol.ndim == 3 and dm_promol.shape[0] == 2:
                dm_promol = np.asarray(dm_promol, dtype=np.float64)
                dm_promol_alpha = dm_promol[0]
                dm_promol_beta = dm_promol[1]
                dm_promol_total = dm_promol_alpha + dm_promol_beta
            else:
                raise ValueError(f"Unexpected promolecule density format: {type(dm_promol)}, shape: {dm_promol.shape if hasattr(dm_promol, 'shape') else 'N/A'}")

            promol_a_e = np.einsum('ij,ji->', dm_promol_alpha, s1e).item()
            promol_b_e = np.einsum('ij,ji->', dm_promol_beta, s1e).item()
            promol_e = promol_a_e + promol_b_e

            print(f"  System: UKS (open-shell)")
            print(f"  Promolecule alpha electrons (Tr[DM*S]): {promol_a_e:.6f}")
            print(f"  Promolecule beta electrons (Tr[DM*S]): {promol_b_e:.6f}")
            print(f"  Promolecule total electrons (Tr[DM*S]): {promol_e:.6f}")

            dm_promol = dm_promol_total

        # Calculate deformation density
        deformation_density = dm_total - dm_promol

        # Verify shapes match
        if deformation_density.shape != (nao, nao):
            raise ValueError(f"Deformation density shape {deformation_density.shape} doesn't match AO basis {nao}x{nao}")

        # Calculate deformation integral (should be close to 0)
        deformation_integral = np.einsum('ij,ji->', deformation_density, s1e).item()
        max_deformation = np.max(np.abs(deformation_density))

        print(f"  Deformation integral (Tr[ΔDM*S]): {deformation_integral:.6f} (should be ~0)")
        print(f"  Max |deformation|: {max_deformation:.6f}")

        # Save cube files
        scf_density_file = os.path.join(OUTPUT_DIR, 'scf_density.cube')
        promol_density_file = os.path.join(OUTPUT_DIR, 'promolecule_density.cube')
        deform_density_file = os.path.join(OUTPUT_DIR, 'deformation_density.cube')

        cubegen.density(mol, scf_density_file, dm_total, nx=nx_g, ny=ny_g, nz=nz_g, margin=margin_bohr)
        print(f"  ✓ SCF density: {scf_density_file}")

        cubegen.density(mol, promol_density_file, dm_promol, nx=nx_g, ny=ny_g, nz=nz_g, margin=margin_bohr)
        print(f"  ✓ Promolecule density: {promol_density_file}")

        cubegen.density(mol, deform_density_file, deformation_density, nx=nx_g, ny=ny_g, nz=nz_g, margin=margin_bohr)
        print(f"  ✓ Deformation density: {deform_density_file}")

        print(f"  Interpretation:")
        print(f"    - Positive regions (red): electron accumulation (bonding)")
        print(f"    - Negative regions (blue): electron depletion (atomic cores)")
    
    print("="*70)

# ============================================================================
# 3. TDDFT CALCULATION
# ============================================================================

if not ENABLE_TDDFT or NUM_EXCITED_STATES <= 0:
    # Skip TDDFT - DFT only mode
    print("\n" + "="*70)
    print("TDDFT SKIPPED")
    print("="*70)
    if not ENABLE_TDDFT:
        print("Reason: ENABLE_TDDFT = False")
    else:
        print("Reason: NUM_EXCITED_STATES = 0")
    print("\n✓ Calculation completed (DFT only)")
    print(f"✓ All files saved to: {OUTPUT_DIR}/\n")
    import sys
    sys.exit(0)

print("\n" + "="*70)
print("TDDFT CALCULATION")
print("="*70)
print(f"Note: TDDFT inherits XC functional ({XC_FUNCTIONAL}) and basis set ({BASIS_SET}) from ground state")
method_name = 'TDA' if USE_TDA else 'TDDFT'
spin_type = 'RKS' if actual_spin == 1 else 'UKS'
print(f"TDDFT method: {method_name} ({spin_type}-based)")

# Select TDA or full TDDFT
# TDDFT automatically uses correct method based on mf (RKS or UKS)
if USE_TDA:
    td = tddft.TDA(mf)  # Faster, ~95% accuracy
    print("Using TDA (Tamm-Dancoff Approximation) - faster calculation")
else:
    td = tddft.TDDFT(mf)  # More accurate
    print("Using full TDDFT - more accurate but slower")

# Set verbosity for TDDFT - use higher verbose to see iteration progress
if VERBOSE_LEVEL >= 3:
    td.verbose = 6  # Debug: show everything including each Davidson iteration
elif VERBOSE_LEVEL >= 2:
    td.verbose = 5  # Detailed: show Davidson iterations
elif VERBOSE_LEVEL >= 1:
    td.verbose = 4  # Normal: show key info
else:
    td.verbose = 3  # Minimal

# Set convergence parameters for TDDFT eigenvalue solver
td.max_cycle = 100  # Maximum Davidson iterations (default is 100)
td.conv_tol = 1e-5  # Convergence tolerance (default is 1e-5)

td.nstates = NUM_EXCITED_STATES
print(f"Calculating {NUM_EXCITED_STATES} excited states...")
print(f"  Max Davidson cycles: {td.max_cycle}")
print(f"  Convergence tolerance: {td.conv_tol}")
print(f"  Verbose level: {td.verbose}")
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

# Print excitation energies with precision based on verbose level
# These are ABSORPTION energies (vertical excitation from ground state geometry)
print("\n" + "="*70)
print("ABSORPTION ENERGIES (Vertical Excitation)")
print("="*70)
print("Ground state geometry → Excited state (Franck-Condon)")
absorption_energies = {}
for i, energy in enumerate(td.e):
    absorption_energies[i] = energy * 27.211  # Store in eV
    if VERBOSE_LEVEL >= 2:
        print(f"State S{i+1}: {energy:.8f} a.u. = {energy*27.211:.6f} eV")
    else:
        print(f"State S{i+1}: {energy:.6f} a.u. = {energy*27.211:.3f} eV")
print("="*70)

# Store ground state total energy for emission calculation
ground_state_energy_gs_geom = mf.e_tot

# ============================================================================
# 3B. EMISSION CALCULATION (Excited State Geometry Optimization)
# ============================================================================

emission_energies = {}
emission_mol = None  # Will store the optimized excited state geometry
td_emission = None
mf_emission = None

if ENABLE_EMISSION and ENABLE_TDDFT:
    print("\n" + "="*70)
    print("EMISSION CALCULATION")
    print("="*70)
    print(f"Optimizing geometry of excited state S{EMISSION_STATE+1}...")
    print(f"This calculates the emission (fluorescence) energy")
    print("-"*70)
    
    # Validate emission state
    if EMISSION_STATE >= NUM_EXCITED_STATES:
        print(f"ERROR: EMISSION_STATE={EMISSION_STATE} but only {NUM_EXCITED_STATES} states calculated")
        print(f"Setting EMISSION_STATE to {NUM_EXCITED_STATES-1}")
        EMISSION_STATE = NUM_EXCITED_STATES - 1
    
    # Import geometry optimization
    from pyscf.geomopt.geometric_solver import optimize as geom_optimize
    
    # Create TDDFT gradient scanner for excited state optimization
    print(f"Setting up excited state gradient for S{EMISSION_STATE+1}...")
    
    # Create gradient scanner for the specified excited state
    # state parameter is 1-indexed in PySCF gradient
    excited_grad = td.nuc_grad_method().as_scanner(state=EMISSION_STATE+1)
    
    print(f"Optimizing S{EMISSION_STATE+1} geometry...")
    print(f"  Max steps: {EMISSION_OPT_MAX_STEPS}")
    print(f"  Convergence: {EMISSION_OPT_CONV}")
    
    # Get convergence parameters
    conv_params = get_opt_conv_params(EMISSION_OPT_CONV)
    
    # Optimize excited state geometry
    emission_mol = geom_optimize(excited_grad, maxsteps=EMISSION_OPT_MAX_STEPS, **conv_params)
    
    # Check optimization convergence and warn if not converged
    opt_converged = getattr(emission_mol, 'converged', None)
    if opt_converged is False:
        print("\n" + "!"*70)
        print("⚠️  WARNING: Excited state geometry optimization did NOT converge!")
        print("   Emission energies and Stokes shifts may be inaccurate!")
        print("   Consider increasing EMISSION_OPT_MAX_STEPS or relaxing convergence criteria.")
        print("!"*70 + "\n")
    elif opt_converged is None:
        print("  Note: Could not verify optimization convergence status.")
    
    # Save optimized excited state geometry
    excited_xyz = os.path.join(OUTPUT_DIR, f'excited_state_S{EMISSION_STATE+1}_geometry.xyz')
    emission_mol.tofile(excited_xyz, format='xyz')
    print(f"✓ Excited state geometry saved to: {excited_xyz}")
    
    # Now calculate emission energy at the optimized excited state geometry
    print("\n" + "-"*70)
    print("Calculating emission energy at excited state geometry...")
    print("-"*70)
    
    # Build DFT at excited state geometry
    # IMPORTANT: Use clean_xc and disp_suffix for consistency with ground state
    if actual_spin == 1:
        mf_emission = dft.RKS(emission_mol, xc=clean_xc)
    else:
        mf_emission = dft.UKS(emission_mol, xc=clean_xc)
    # Apply dispersion correction if used in ground state (CRITICAL for consistency)
    if disp_suffix:
        mf_emission.disp = disp_suffix
    # Use same grid level as ground state mf for consistency
    mf_emission.grids.level = mf.grids.level
    mf_emission.verbose = 0
    mf_emission.kernel()
    
    # Calculate TDDFT at excited state geometry
    if USE_TDA:
        td_emission = tddft.TDA(mf_emission)
    else:
        td_emission = tddft.TDDFT(mf_emission)
    td_emission.nstates = NUM_EXCITED_STATES
    td_emission.verbose = 0
    td_emission.kernel()
    
    # Emission energy is the excitation energy at the excited state geometry
    # This represents the vertical de-excitation (fluorescence)
    print("\n" + "="*70)
    print("EMISSION ENERGIES (Vertical De-excitation)")
    print("="*70)
    print("Excited state geometry → Ground state (Fluorescence)")
    for i, energy in enumerate(td_emission.e):
        emission_energies[i] = energy * 27.211  # Store in eV
        if VERBOSE_LEVEL >= 2:
            print(f"State S{i+1}: {energy:.8f} a.u. = {energy*27.211:.6f} eV")
        else:
            print(f"State S{i+1}: {energy:.6f} a.u. = {energy*27.211:.3f} eV")
    print("="*70)
    
    # Calculate Stokes shift
    print("\n" + "="*70)
    print("STOKES SHIFT (Absorption - Emission)")
    print("="*70)
    print(f"⚠️  NOTE: Only S{EMISSION_STATE+1} Stokes shift is physically meaningful!")
    print(f"   Geometry was optimized for S{EMISSION_STATE+1} only.")
    print(f"   Other states' Stokes shifts are calculated at S{EMISSION_STATE+1} geometry.")
    print("-"*70)
    for i in range(min(len(absorption_energies), len(emission_energies))):
        stokes = absorption_energies[i] - emission_energies[i]
        validity = "✓ (valid)" if i == EMISSION_STATE else "(at S{} geom)".format(EMISSION_STATE+1)
        print(f"State S{i+1}: {stokes:.6f} eV = {stokes*8065.544:.1f} cm⁻¹  {validity}")
    print("="*70)
    
    # Summary table
    print("\n" + "="*70)
    print("ABSORPTION vs EMISSION SUMMARY")
    print("="*70)
    print(f"{'State':<8} {'Absorption (eV)':<18} {'Emission (eV)':<18} {'Stokes (eV)':<15} {'Stokes (cm⁻¹)':<15}")
    print("-"*70)
    for i in range(min(len(absorption_energies), len(emission_energies))):
        abs_e = absorption_energies[i]
        emi_e = emission_energies[i]
        stokes = abs_e - emi_e
        stokes_cm = stokes * 8065.544
        print(f"S{i+1:<7} {abs_e:<18.6f} {emi_e:<18.6f} {stokes:<15.6f} {stokes_cm:<15.1f}")
    print("="*70)

elif ENABLE_EMISSION and not ENABLE_TDDFT:
    print("\nNote: ENABLE_EMISSION=True but ENABLE_TDDFT=False. Skipping emission calculation.")

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

# Calculate dipole integrals (int1e_r is the modern PySCF interface)
dip_ints = mol.intor('int1e_r', comp=3)  # x, y, z components

def calculate_transition_dipole(td, state_id):
    """
    Calculate transition dipole moment between ground state and excited state.
    Handles both RKS and UKS systems.
    
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
    
    # Check if TDA (Y=0 as integer, not array)
    is_tda = not hasattr(Y, 'shape')
    
    if is_uks:
        # UKS: mo_coeff and mo_occ are tuples (alpha, beta)
        # For UKS, X and Y are also tuples: (Xa, Xb), (Ya, Yb)
        mo_coeff_a, mo_coeff_b = mo_coeff
        mo_occ_a, mo_occ_b = mo_occ
        Xa, Xb = X
        if is_tda:
            Ya = Yb = 0
        else:
            Ya, Yb = Y
        
        # Alpha spin (for TDA: Y=0, so X+Y = X)
        nocc_a = Xa.shape[0]
        nvir_a = Xa.shape[1]
        nmo_a = mo_coeff_a.shape[1]
        
        amp_a = Xa if is_tda else (Xa + Ya)
        t_dm1_mo_a = np.zeros((nmo_a, nmo_a))
        t_dm1_mo_a[:nocc_a, nocc_a:] = amp_a.reshape(nocc_a, nvir_a)
        t_dm1_ao_a = reduce(np.dot, (mo_coeff_a, t_dm1_mo_a, mo_coeff_a.T))
        
        # Beta spin
        nocc_b = Xb.shape[0]
        nvir_b = Xb.shape[1]
        nmo_b = mo_coeff_b.shape[1]
        
        amp_b = Xb if is_tda else (Xb + Yb)
        t_dm1_mo_b = np.zeros((nmo_b, nmo_b))
        t_dm1_mo_b[:nocc_b, nocc_b:] = amp_b.reshape(nocc_b, nvir_b)
        t_dm1_ao_b = reduce(np.dot, (mo_coeff_b, t_dm1_mo_b, mo_coeff_b.T))
        
        # Total transition density (alpha + beta)
        t_dm1_ao = t_dm1_ao_a + t_dm1_ao_b
        
    else:
        # RKS: mo_coeff and mo_occ are arrays
        orbo = mo_coeff[:, mo_occ > 0]
        orbv = mo_coeff[:, mo_occ == 0]
        nocc = orbo.shape[1]
        nvir = orbv.shape[1]
        
        # Transition density in MO basis (for TDA: Y=0, so X+Y = X)
        amp = X if is_tda else (X + Y)
        t_dm1_mo = np.zeros((mo_coeff.shape[1], mo_coeff.shape[1]))
        t_dm1_mo[:nocc, nocc:] = amp.reshape(nocc, nvir)
        
        # Transform to AO basis
        t_dm1_ao = reduce(np.dot, (mo_coeff, t_dm1_mo, mo_coeff.T))
    
    # Calculate transition dipole: μ = Tr(μ_op * T)
    tdm = np.einsum('xij,ji->x', dip_ints, t_dm1_ao)
    
    return tdm

def calculate_transition_dipole_with_ints(td_obj, state_id, dip_ints_use):
    tdm_ao = calculate_transition_dipole(td_obj, state_id)
    return tdm_ao, float(np.linalg.norm(tdm_ao))


# Calculate and print transition dipoles for all states
print("\n--- ABSORPTION Transition Dipole Moments (at ground state geometry) ---")
print("These TDMs correspond to S₀ → Sₙ transitions (light absorption)")
print(f"\n{'State':<8} {'μ_x':<12} {'μ_y':<12} {'μ_z':<12} {'|μ|':<12} {'f':<12}")
print("-" * 70)

for i in range(td.nstates):
    tdm, tdm_magnitude = calculate_transition_dipole_with_ints(td, i, dip_ints)

    # Oscillator strength: f = (2/3) * ω * |μ|^2
    # where ω is excitation energy in a.u.
    omega = td.e[i]
    osc_strength = (2.0/3.0) * omega * tdm_magnitude**2

    print(f"{i+1:<8} {tdm[0]:>11.6f} {tdm[1]:>11.6f} {tdm[2]:>11.6f} "
          f"{tdm_magnitude:>11.6f} {osc_strength:>11.6f}")

# EMISSION TDMs: at excited state geometry (if emission calculation was done)
if td_emission is not None and emission_mol is not None:
    print("\n--- EMISSION Transition Dipole Moments (at excited state geometry) ---")
    print("These TDMs correspond to Sₙ → S₀ transitions (fluorescence emission)")
    print("Note: Calculated from TDDFT at the optimized excited state geometry")

    emission_charges = emission_mol.atom_charges()
    emission_coords = emission_mol.atom_coords()
    emission_nuc_center = np.einsum('z,zx->x', emission_charges, emission_coords) / emission_charges.sum()
    emission_mol.set_common_orig_(emission_nuc_center)
    dip_ints_emission = emission_mol.intor('int1e_r', comp=3)

    print(f"\n{'State':<8} {'μ_x':<12} {'μ_y':<12} {'μ_z':<12} {'|μ|':<12} {'f':<12}")
    print("-" * 70)

    for i in range(td_emission.nstates):
        tdm_em, tdm_em_magnitude = calculate_transition_dipole_with_ints(td_emission, i, dip_ints_emission)
        omega_em = td_emission.e[i]
        osc_strength_em = (2.0/3.0) * omega_em * tdm_em_magnitude**2
        print(f"{i+1:<8} {tdm_em[0]:>11.6f} {tdm_em[1]:>11.6f} {tdm_em[2]:>11.6f} "
              f"{tdm_em_magnitude:>11.6f} {osc_strength_em:>11.6f}")

print("="*70)

# ============================================================================
# 5. TRANSITION DENSITY MATRICES
# ============================================================================

def calculate_transition_density_matrix(td, state_id):
    """
    Calculate transition density matrix between ground and excited state.
    This is the proper transition density for visualization.
    Handles both RKS and UKS systems.
    
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
    
    # Check if TDA (Y=0 as integer, not array)
    is_tda = not hasattr(Y, 'shape')
    
    if is_uks:
        mo_coeff_a, mo_coeff_b = mo_coeff
        mo_occ_a, mo_occ_b = mo_occ
        Xa, Xb = X
        if is_tda:
            Ya = Yb = 0
        else:
            Ya, Yb = Y
        
        nocc_a = Xa.shape[0]
        nvir_a = Xa.shape[1]
        nocc_b = Xb.shape[0]
        nvir_b = Xb.shape[1]
        nmo_a = mo_coeff_a.shape[1]
        nmo_b = mo_coeff_b.shape[1]
        
        # Alpha spin (for TDA: Y=0, so X+Y = X)
        amp_a = Xa if is_tda else (Xa + Ya)
        t_dm1_mo_a = np.zeros((nmo_a, nmo_a))
        t_dm1_mo_a[:nocc_a, nocc_a:] = amp_a
        t_dm1_mo_a[nocc_a:, :nocc_a] = amp_a.T
        t_dm1_ao_a = reduce(np.dot, (mo_coeff_a, t_dm1_mo_a, mo_coeff_a.T))
        
        # Beta spin
        amp_b = Xb if is_tda else (Xb + Yb)
        t_dm1_mo_b = np.zeros((nmo_b, nmo_b))
        t_dm1_mo_b[:nocc_b, nocc_b:] = amp_b
        t_dm1_mo_b[nocc_b:, :nocc_b] = amp_b.T
        t_dm1_ao_b = reduce(np.dot, (mo_coeff_b, t_dm1_mo_b, mo_coeff_b.T))
        
        t_dm1_ao = t_dm1_ao_a + t_dm1_ao_b
    else:
        orbo = mo_coeff[:, mo_occ > 0]
        orbv = mo_coeff[:, mo_occ == 0]
        nocc = orbo.shape[1]
        nvir = orbv.shape[1]
        
        # Transition density in MO basis (for TDA: Y=0, so X+Y = X)
        amp = X.reshape(nocc, nvir) if is_tda else (X + Y).reshape(nocc, nvir)
        t_dm1_mo = np.zeros((mo_coeff.shape[1], mo_coeff.shape[1]))
        t_dm1_mo[:nocc, nocc:] = amp
        t_dm1_mo[nocc:, :nocc] = amp.T
        t_dm1_ao = reduce(np.dot, (mo_coeff, t_dm1_mo, mo_coeff.T))
    
    return t_dm1_ao

def calculate_excited_state_density(td, state_id):
    """
    Calculate excited state density matrix.
    Based on PySCF examples/tddft/22-density.py
    Handles both RKS and UKS systems.
    
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
    
    # Check if TDA (Y=0 as integer, not array)
    is_tda = not hasattr(Y, 'shape')
    
    if is_uks:
        # UKS case - proper excited state density
        # ρ_excited = ρ_ground + Δρ_alpha + Δρ_beta
        # where Δρ comes from the TDDFT response (X,Y amplitudes)
        Xa, Xb = X

        # For UKS, Y can be (Ya,Yb) or 0 (TDA)
        is_tda = not isinstance(Y, tuple) or (isinstance(Y, tuple) and not hasattr(Y[0], 'shape'))
        if is_tda:
            Ya = Yb = 0
        else:
            Ya, Yb = Y

        # mo_coeff/mo_occ can be tuples (alpha,beta)
        mo_coeff_a, mo_coeff_b = mo_coeff
        mo_occ_a, mo_occ_b = mo_occ

        dm0 = mf.make_rdm1()
        if isinstance(dm0, tuple):
            dm0_a, dm0_b = dm0
        else:
            dm0 = np.asarray(dm0)
            if dm0.ndim == 3 and dm0.shape[0] == 2:
                dm0_a, dm0_b = dm0[0], dm0[1]
            else:
                raise ValueError(f"Unexpected UKS ground-state density matrix format: {type(dm0)}, shape={getattr(dm0,'shape',None)}")

        occ_a = mo_occ_a > 0
        vir_a = mo_occ_a == 0
        orbo_a = mo_coeff_a[:, occ_a]
        orbv_a = mo_coeff_a[:, vir_a]
        if Xa.shape != (orbo_a.shape[1], orbv_a.shape[1]):
            raise ValueError(f"Xa shape {Xa.shape} does not match (nocc,nvir)=({orbo_a.shape[1]},{orbv_a.shape[1]})")

        dm_oo_a = -np.einsum('ia,ka->ik', Xa.conj(), Xa)
        dm_vv_a = np.einsum('ia,ic->ac', Xa, Xa.conj())
        if not is_tda:
            dm_oo_a -= np.einsum('ia,ka->ik', Ya.conj(), Ya)
            dm_vv_a += np.einsum('ia,ic->ac', Ya, Ya.conj())

        delta_ao_a = (
            np.einsum('pi,ij,qj->pq', orbo_a, dm_oo_a, orbo_a.conj())
            + np.einsum('pa,ab,qb->pq', orbv_a, dm_vv_a, orbv_a.conj())
        )
        dm_ao_a = np.asarray(dm0_a, dtype=np.float64) + np.asarray(delta_ao_a.real, dtype=np.float64)

        occ_b = mo_occ_b > 0
        vir_b = mo_occ_b == 0
        orbo_b = mo_coeff_b[:, occ_b]
        orbv_b = mo_coeff_b[:, vir_b]
        if Xb.shape != (orbo_b.shape[1], orbv_b.shape[1]):
            raise ValueError(f"Xb shape {Xb.shape} does not match (nocc,nvir)=({orbo_b.shape[1]},{orbv_b.shape[1]})")

        dm_oo_b = -np.einsum('ia,ka->ik', Xb.conj(), Xb)
        dm_vv_b = np.einsum('ia,ic->ac', Xb, Xb.conj())
        if not is_tda:
            dm_oo_b -= np.einsum('ia,ka->ik', Yb.conj(), Yb)
            dm_vv_b += np.einsum('ia,ic->ac', Yb, Yb.conj())

        delta_ao_b = (
            np.einsum('pi,ij,qj->pq', orbo_b, dm_oo_b, orbo_b.conj())
            + np.einsum('pa,ab,qb->pq', orbv_b, dm_vv_b, orbv_b.conj())
        )
        dm_ao_b = np.asarray(dm0_b, dtype=np.float64) + np.asarray(delta_ao_b.real, dtype=np.float64)

        dm_excited = dm_ao_a + dm_ao_b
        return np.asarray(dm_excited.real, dtype=np.float64)
    else:
        # RKS case
        # For TDDFT, need to consider both X and Y
        # Density matrix changes in MO basis
        nocc = X.shape[0]
        
        # Occupied-occupied and virtual-virtual blocks
        dm_oo = -np.einsum('ia,ka->ik', X.conj(), X)
        dm_vv = np.einsum('ia,ic->ac', X, X.conj())
        
        # Add Y contribution only for full TDDFT (not TDA)
        if not is_tda:
            dm_oo -= np.einsum('ia,ka->ik', Y.conj(), Y)
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

def get_orbital_labels_single_spin(mo_occ, mo_coeff, spin_label='', mo_energy=None):
    """
    Get orbital labels for a single spin channel.
    Based on PySCF's official implementation.
    
    Args:
        mo_occ: occupation numbers for this spin
        mo_coeff: MO coefficients for this spin  
        spin_label: 'α' or 'β' for labeling
        mo_energy: orbital energies (optional, for enhanced output)
    
    Returns:
        labels: list of orbital labels
        homo_idx: index of HOMO (0-indexed)
        nocc: number of occupied orbitals
        energies_eV: orbital energies in eV (or None if not provided)
    """
    nmo = mo_coeff.shape[1]
    homo_idx = np.where(mo_occ > 0)[0][-1]
    nocc = homo_idx + 1
    
    # Convert energies to eV if provided
    energies_eV = None
    if mo_energy is not None:
        energies_eV = np.asarray(mo_energy) * 27.211  # Hartree to eV
    
    labels = []
    for i in range(nmo):
        if i <= homo_idx:
            offset = homo_idx - i
            if offset == 0:
                label = f'HOMO({spin_label})' if spin_label else 'HOMO'
            else:
                label = f'HOMO-{offset}({spin_label})' if spin_label else f'HOMO-{offset}'
        else:
            offset = i - homo_idx - 1
            if offset == 0:
                label = f'LUMO({spin_label})' if spin_label else 'LUMO'
            else:
                label = f'LUMO+{offset}({spin_label})' if spin_label else f'LUMO+{offset}'
        labels.append(label)
    
    return labels, homo_idx, nocc, energies_eV


def get_orbital_labels(mf):
    """
    Get orbital labels for cube file naming.
    Works for both RKS and UKS systems.
    For UKS, uses alpha spin labels (without spin designation for file naming).
    
    Returns:
        labels: list of orbital labels (HOMO, LUMO, HOMO-1, etc.)
        homo_idx: index of HOMO
    """
    mo_occ = mf.mo_occ
    mo_coeff = mf.mo_coeff
    
    # Handle UKS (use alpha spin)
    if isinstance(mo_occ, tuple):
        mo_occ = mo_occ[0]
    if isinstance(mo_coeff, tuple):
        mo_coeff = mo_coeff[0]
    
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


def analyze_single_spin_contribution(X, Y, labels, spin_label, threshold=0.01, energies_eV=None):
    """
    Analyze contributions from a single spin channel.
    Based on PySCF's official UHF TDDFT implementation.
    
    Args:
        X: TDDFT X amplitudes (nocc x nvir)
        Y: TDDFT Y amplitudes (nocc x nvir) or 0 for TDA
        labels: orbital labels for this spin
        spin_label: 'α' or 'β'
        threshold: minimum weight to include
        energies_eV: orbital energies in eV (optional)
    
    Returns:
        list of (occ_idx, vir_idx, weight, label, spin, occ_energy, vir_energy)
    """
    # Check if TDA
    is_tda = not hasattr(Y, 'shape')
    
    nocc, nvir = X.shape
    
    # Calculate amplitudes: X+Y for full TDDFT, X for TDA
    amplitudes = X if is_tda else (X + Y)
    
    # Weights are squared amplitudes
    weights = amplitudes ** 2
    
    contributions = []
    for i in range(nocc):
        for a in range(nvir):
            weight = weights[i, a]
            if weight > threshold:
                occ_idx = i  # 0-indexed occupied orbital
                vir_idx = nocc + a  # 0-indexed virtual orbital (absolute index)
                occ_label = labels[occ_idx]
                vir_label = labels[vir_idx]
                transition_label = f"{occ_label} → {vir_label}"
                
                # Get orbital energies if available
                occ_energy = energies_eV[occ_idx] if energies_eV is not None else None
                vir_energy = energies_eV[vir_idx] if energies_eV is not None else None
                
                contributions.append((occ_idx, vir_idx, weight, transition_label, spin_label, occ_energy, vir_energy))
    
    return contributions


def analyze_transition_contributions(td, state_id, mf, threshold=0.01, top_n=10):
    """
    Analyze orbital pair contributions to a specific excited state.
    Properly handles both RKS and UKS systems following PySCF's official implementation.
    
    Args:
        td: TDDFT object
        state_id: excited state index (0-indexed)
        mf: mean-field object
        threshold: minimum weight to include (default 0.01 = 1%)
        top_n: number of top contributions to return
    
    Returns:
        contributions: list of (occ_idx, vir_idx, weight, label, spin, occ_energy, vir_energy) sorted by weight
        total_weight: sum of all weights
    """
    X, Y = td.xy[state_id]
    
    mo_occ = mf.mo_occ
    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy
    
    # Check if UKS (unrestricted) - for UKS, X is a tuple of (Xa, Xb)
    # This is more reliable than checking mo_occ which may vary by implementation
    is_uks = isinstance(X, tuple) and len(X) == 2 and hasattr(X[0], 'shape')
    
    # For TDA, Y is 0 (scalar) or tuple of zeros
    if is_uks:
        # Check if TDA by looking at Y structure
        is_tda = not isinstance(Y, tuple) or (isinstance(Y, tuple) and not hasattr(Y[0], 'shape'))
    else:
        is_tda = not hasattr(Y, 'shape')
    
    all_contributions = []
    
    if is_uks:
        # UKS: analyze both spins
        # Handle mo_occ/mo_coeff/mo_energy which can be tuple OR 2D array
        if isinstance(mo_occ, tuple):
            mo_occ_a, mo_occ_b = mo_occ
        else:
            mo_occ_a, mo_occ_b = mo_occ[0], mo_occ[1]
        
        if isinstance(mo_coeff, tuple):
            mo_coeff_a, mo_coeff_b = mo_coeff
        else:
            mo_coeff_a, mo_coeff_b = mo_coeff[0], mo_coeff[1]
        
        if isinstance(mo_energy, tuple):
            mo_energy_a, mo_energy_b = mo_energy
        else:
            mo_energy_a, mo_energy_b = mo_energy[0], mo_energy[1]
        
        Xa, Xb = X
        Ya, Yb = (0, 0) if is_tda else Y
        
        # Get labels for both spins (with energies)
        labels_a, homo_a, nocc_a, energies_a = get_orbital_labels_single_spin(mo_occ_a, mo_coeff_a, 'α', mo_energy_a)
        labels_b, homo_b, nocc_b, energies_b = get_orbital_labels_single_spin(mo_occ_b, mo_coeff_b, 'β', mo_energy_b)
        
        # Analyze alpha spin
        contrib_a = analyze_single_spin_contribution(Xa, Ya, labels_a, 'α', threshold, energies_a)
        all_contributions.extend(contrib_a)
        
        # Analyze beta spin
        contrib_b = analyze_single_spin_contribution(Xb, Yb, labels_b, 'β', threshold, energies_b)
        all_contributions.extend(contrib_b)
        
    else:
        # RKS: single spin channel
        labels, homo_idx, nocc, energies = get_orbital_labels_single_spin(mo_occ, mo_coeff, '', mo_energy)
        contrib = analyze_single_spin_contribution(X, Y, labels, '', threshold, energies)
        all_contributions.extend(contrib)
    
    # Sort by weight (descending)
    all_contributions.sort(key=lambda x: x[2], reverse=True)
    
    # Normalize weights to sum to 1
    total_weight = sum(c[2] for c in all_contributions)
    if total_weight > 0:
        all_contributions = [(occ, vir, w/total_weight, label, spin, occ_e, vir_e) 
                             for occ, vir, w, label, spin, occ_e, vir_e in all_contributions]
    
    return all_contributions[:top_n], total_weight

def calculate_pair_transition_density(mf, occ_idx, vir_idx, spin='alpha'):
    """
    Calculate transition density matrix for a single orbital pair i→a.
    Handles both RKS and UKS (uses alpha spin for UKS).
    
    T_μν = C_μ^i × C_ν^a + C_μ^a × C_ν^i
    
    Returns NumPy array in AO basis.
    """
    mo_coeff = mf.mo_coeff
    
    # Ensure we have NumPy array
    mo_coeff = np.asarray(mo_coeff)
    
    is_uks = False
    if isinstance(mo_coeff, tuple):
        is_uks = True
        mo_coeff_alpha = np.asarray(mo_coeff[0])
        mo_coeff_beta = np.asarray(mo_coeff[1])
    elif mo_coeff.ndim == 3 and mo_coeff.shape[0] == 2:
        is_uks = True
        mo_coeff_alpha = np.asarray(mo_coeff[0])
        mo_coeff_beta = np.asarray(mo_coeff[1])
    else:
        mo_coeff_use = mo_coeff

    if is_uks:
        if spin.lower() == 'beta':
            mo_coeff_use = mo_coeff_beta
        else:
            mo_coeff_use = mo_coeff_alpha

    mo_coeff_use = np.asarray(mo_coeff_use)
    if mo_coeff_use.ndim != 2:
        raise ValueError(f"mo_coeff must be 2D, got shape {mo_coeff_use.shape}, ndim={mo_coeff_use.ndim}")
    
    # Check if indices are valid
    nao, nmo = mo_coeff_use.shape
    if occ_idx >= nmo or vir_idx >= nmo:
        raise IndexError(f"Orbital indices out of range: occ_idx={occ_idx}, vir_idx={vir_idx}, but only {nmo} MOs available (nao={nao})")
    
    # Extract specific orbitals
    occ_mo = mo_coeff_use[:, occ_idx]
    vir_mo = mo_coeff_use[:, vir_idx]
    
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
            
            # Print contribution table with orbital indices and energies
            print(f"\n{'='*100}")
            print(f"STATE {state_id+1}: {excitation_energy:.4f} eV")
            print(f"{'='*100}")
            print(f"{'Rank':<5} {'Idx':<8} {'Transition':<28} {'Occ E(eV)':<12} {'Vir E(eV)':<12} {'ΔE(eV)':<10} {'Weight':<10} {'%':<8}")
            print(f"{'-'*100}")
            
            cumulative = 0.0
            for rank, (occ_idx, vir_idx, weight, label, spin, occ_e, vir_e) in enumerate(contributions, 1):
                cumulative += weight
                idx_str = f"{occ_idx}→{vir_idx}"
                delta_e = (vir_e - occ_e) if (occ_e is not None and vir_e is not None) else 0.0
                occ_e_str = f"{occ_e:.3f}" if occ_e is not None else "N/A"
                vir_e_str = f"{vir_e:.3f}" if vir_e is not None else "N/A"
                delta_e_str = f"{delta_e:.3f}" if (occ_e is not None and vir_e is not None) else "N/A"
                print(f"{rank:<5} {idx_str:<8} {label:<28} {occ_e_str:<12} {vir_e_str:<12} {delta_e_str:<10} {weight:<10.4f} {weight*100:<7.2f}%")
            
            print(f"{'-'*100}")
            print(f"Total weight analyzed: {total_weight:.6f}")
            
            # Note about orbital energies
            print(f"\nNote: Idx = orbital indices (0-indexed), E = orbital energy in eV")
            print(f"      ΔE = single-particle energy gap (Kohn-Sham), not excitation energy")
        
        # Save contribution tables to file (with orbital indices and energies)
        table_file = os.path.join(OUTPUT_DIR, 'contribution_tables.txt')
        with open(table_file, 'w') as f:
            f.write("="*100 + "\n")
            f.write("ORBITAL PAIR CONTRIBUTIONS TO EXCITED STATES\n")
            f.write("="*100 + "\n\n")
            
            for state_id in valid_contrib_states:
                contributions, total_weight = all_contributions[state_id]
                excitation_energy = td.e[state_id] * 27.211
                
                f.write(f"\n{'='*100}\n")
                f.write(f"STATE {state_id+1}: {excitation_energy:.4f} eV\n")
                f.write(f"{'='*100}\n")
                f.write(f"{'Rank':<5} {'Idx':<8} {'Transition':<28} {'Occ E(eV)':<12} {'Vir E(eV)':<12} {'ΔE(eV)':<10} {'Weight':<10} {'%':<8}\n")
                f.write(f"{'-'*100}\n")
                
                cumulative = 0.0
                for rank, (occ_idx, vir_idx, weight, label, spin, occ_e, vir_e) in enumerate(contributions, 1):
                    cumulative += weight
                    idx_str = f"{occ_idx}→{vir_idx}"
                    delta_e = (vir_e - occ_e) if (occ_e is not None and vir_e is not None) else 0.0
                    occ_e_str = f"{occ_e:.3f}" if occ_e is not None else "N/A"
                    vir_e_str = f"{vir_e:.3f}" if vir_e is not None else "N/A"
                    delta_e_str = f"{delta_e:.3f}" if (occ_e is not None and vir_e is not None) else "N/A"
                    f.write(f"{rank:<5} {idx_str:<8} {label:<28} {occ_e_str:<12} {vir_e_str:<12} {delta_e_str:<10} {weight:<10.4f} {weight*100:<7.2f}%\n")
                
                f.write(f"{'-'*100}\n")
                f.write(f"Total weight analyzed: {total_weight:.6f}\n")
                f.write(f"Note: Idx = orbital indices (0-indexed), E = orbital energy in eV\n")
                f.write(f"      ΔE = single-particle gap (Kohn-Sham), Contributions = (amplitude)²\n")
                f.write(f"{'='*100}\n\n")
        
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
    
    homo_idx = np.where(mo_occ > 0)[0][-1]
    lumo_idx = np.where(mo_occ == 0)[0][0]
    
    print(f"HOMO index: {homo_idx}")
    print(f"LUMO index: {lumo_idx}")
    print(f"HOMO energy: {mo_energy[homo_idx]*27.211:.3f} eV")
    print(f"LUMO energy: {mo_energy[lumo_idx]*27.211:.3f} eV")
    print(f"HOMO-LUMO gap: {(mo_energy[lumo_idx] - mo_energy[homo_idx])*27.211:.3f} eV")
    
    # Generate HOMO cube file
    homo_file = os.path.join(OUTPUT_DIR, 'HOMO.cube')
    cubegen.orbital(mol, homo_file, mo_coeff[:, homo_idx], nx=nx, ny=ny, nz=nz, margin=margin_bohr)
    print(f"\n  ✓ HOMO orbital: {homo_file}")
    
    # Generate LUMO cube file
    lumo_file = os.path.join(OUTPUT_DIR, 'LUMO.cube')
    cubegen.orbital(mol, lumo_file, mo_coeff[:, lumo_idx], nx=nx, ny=ny, nz=nz, margin=margin_bohr)
    print(f"  ✓ LUMO orbital: {lumo_file}")
    
    # Generate HOMO-1 and LUMO+1 for additional verification
    if homo_idx > 0:
        homo1_file = os.path.join(OUTPUT_DIR, 'HOMO-1.cube')
        cubegen.orbital(mol, homo1_file, mo_coeff[:, homo_idx-1], nx=nx, ny=ny, nz=nz, margin=margin_bohr)
        print(f"  ✓ HOMO-1 orbital: {homo1_file}")
    
    if lumo_idx < len(mo_occ) - 1:
        lumo1_file = os.path.join(OUTPUT_DIR, 'LUMO+1.cube')
        cubegen.orbital(mol, lumo1_file, mo_coeff[:, lumo_idx+1], nx=nx, ny=ny, nz=nz, margin=margin_bohr)
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
    cubegen.density(mol, homo_lumo_file, t_homo_lumo, nx=nx, ny=ny, nz=nz, margin=margin_bohr)
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
            cubegen.density(mol, filename_trans, dm_trans, nx=nx, ny=ny, nz=nz, margin=margin_bohr)
            print(f"  ✓ Transition density: {filename_trans}")
        
        # 2. Excited state density
        if GENERATE_EXCITED_DENSITY:
            dm_excited = calculate_excited_state_density(td, state_id)
            filename_excited = os.path.join(OUTPUT_DIR, f'excited_state_density_state{state_id+1}.cube')
            cubegen.density(mol, filename_excited, dm_excited, nx=nx, ny=ny, nz=nz, margin=margin_bohr)
            print(f"  ✓ Excited state density: {filename_excited}")
        
        # 3. Density difference
        if GENERATE_DENSITY_DIFFERENCE:
            if not GENERATE_EXCITED_DENSITY:
                dm_excited = calculate_excited_state_density(td, state_id)
            
            # Get ground state density
            dm_ground = mf.make_rdm1()
            
            # For UKS, dm_ground is a tuple (alpha, beta), sum them
            if isinstance(dm_ground, tuple):
                dm_ground = dm_ground[0] + dm_ground[1]
            
            dm_diff = dm_excited - dm_ground
            filename_diff = os.path.join(OUTPUT_DIR, f'density_difference_state{state_id+1}.cube')
            cubegen.density(mol, filename_diff, dm_diff, nx=nx, ny=ny, nz=nz, margin=margin_bohr)
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
            
            # Check if TDA (Y=0 as integer, not array)
            is_tda = not hasattr(Y, 'shape')
            
            # Handle UKS: X and Y are tuples (Xa, Xb), (Ya, Yb)
            if actual_spin > 1:  # UKS
                Xa, Xb = X
                if is_tda:
                    Ya = 0
                else:
                    Ya, Yb = Y
                # Use alpha for verification
                nocc = Xa.shape[0]
                nvir = Xa.shape[1]
                homo_lumo_amplitude = abs(Xa[nocc-1, 0]) if is_tda else abs(Xa[nocc-1, 0] + Ya[nocc-1, 0])
                total_amplitude = np.linalg.norm(Xa) if is_tda else np.linalg.norm(Xa + Ya)
            else:  # RKS
                nocc = X.shape[0]
                nvir = X.shape[1]
                homo_lumo_amplitude = abs(X[nocc-1, 0]) if is_tda else abs(X[nocc-1, 0] + Y[nocc-1, 0])
                total_amplitude = np.linalg.norm(X) if is_tda else np.linalg.norm(X + Y)
            
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
        for rank, (occ_idx, vir_idx, weight, label, spin, occ_e, vir_e) in enumerate(contributions, 1):
            if pair_count >= MAX_PAIRS_PER_STATE:
                break
            
            if weight < PAIR_CONTRIBUTION_CUTOFF:
                print(f"  Skipping {label} (contribution {weight*100:.2f}% < {PAIR_CONTRIBUTION_CUTOFF*100:.0f}%)")
                continue
            
            spin_channel = 'beta' if spin == 'β' else 'alpha'

            # Calculate transition density for this pair
            t_dm_pair = calculate_pair_transition_density(mf, occ_idx, vir_idx, spin=spin_channel)
            
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

            spin_suffix = f"_{spin_channel}" if spin else ""
            
            filename = os.path.join(OUTPUT_DIR, 
                f'transition_pair_state{state_id+1}_{occ_label}_to_{vir_label}{spin_suffix}.cube')

            cubegen.density(mol, filename, t_dm_pair, nx=nx, ny=ny, nz=nz, margin=margin_bohr)
            print(f"  ✓ Rank {rank}: {label} ({weight*100:.2f}%) → {filename}")
            pair_count += 1
    
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
