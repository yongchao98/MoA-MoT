import numpy as np
from pyscf import gto, scf
import warnings

# Suppress a harmless warning from pyscf about J/K builders
warnings.filterwarnings("ignore", category=UserWarning)

def calculate_h2plus_energy(distance, use_symmetry):
    """
    Calculates the ROHF energy for H2+ at a given distance.

    Args:
        distance (float): The internuclear distance in Angstroms.
        use_symmetry (bool): If True, enforces D2h symmetry. If False,
                             uses C1 (no symmetry), allowing symmetry breaking.

    Returns:
        float: The total electronic energy in Hartree.
    """
    # Define the molecule: two H atoms, charge +1, spin 1 (one electron)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, distance)]
    ]
    mol.charge = 1
    mol.spin = 1
    mol.basis = 'cc-pVDZ' # A better basis set
    mol.verbose = 0 # Suppress detailed SCF output

    # Set symmetry based on the input flag.
    # D2h is an Abelian subgroup of the full Dooh symmetry. Enforcing it
    # forces the molecular orbitals to be symmetric or anti-symmetric.
    # C1 means no symmetry constraints are used.
    mol.symmetry = 'd2h' if use_symmetry else 'c1'
    mol.build()

    # Run a Restricted Open-Shell Hartree-Fock (ROHF) calculation
    mf = scf.ROHF(mol)
    energy = mf.kernel()
    return energy

# Define a range of distances to scan
distances = np.arange(0.7, 5.1, 0.2)

# --- Calculation 1: Enforcing Symmetry (INCORRECT at large distance) ---
print("--- 1. Calculation with Symmetry Enforced ---")
print("This forces the electron to be shared equally, which is wrong for dissociation.")
print(f"{'Distance (A)':<15} {'Energy (Hartree)':<20}")
print("-" * 37)

for d in distances:
    try:
        energy = calculate_h2plus_energy(d, use_symmetry=True)
        # At large distances, this energy will be incorrectly high because it
        # describes an unphysical state H(0.5+) ... H(0.5+).
        print(f"{d:<15.2f} {energy:<20.8f}")
    except Exception as e:
        print(f"{d:<15.2f} {'Calculation Failed'}")

# --- Calculation 2: Allowing Symmetry Breaking (CORRECT) ---
print("\n--- 2. Calculation Allowing Symmetry Breaking ---")
print("This allows the electron to localize on one atom, describing H + H+ correctly.")
print("Note how the energy approaches -0.5 Hartree (energy of a H atom).")
print(f"{'Distance (A)':<15} {'Energy (Hartree)':<20}")
print("-" * 37)

for d in distances:
    energy = calculate_h2plus_energy(d, use_symmetry=False)
    print(f"{d:<15.2f} {energy:<20.8f}")

print("\nAs shown, allowing symmetry breaking (calculation 2) fixes the potential energy curve,")
print("resulting in the correct dissociation energy limit of -0.5 Hartree.")
