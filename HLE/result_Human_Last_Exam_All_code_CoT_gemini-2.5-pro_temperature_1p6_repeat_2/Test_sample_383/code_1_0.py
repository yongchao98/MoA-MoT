import numpy as np
try:
    from pyscf import gto, scf
except ImportError:
    print("PySCF is not installed. Please install it using: pip install pyscf")
    exit()

def calculate_h2_plus_energy(distance, allow_symmetry_breaking):
    """
    Calculates the energy of the H2+ cation at a given distance.

    Args:
        distance (float): The internuclear distance in Angstroms.
        allow_symmetry_breaking (bool): If True, symmetry is not enforced,
                                        allowing the wavefunction to localize.

    Returns:
        float: The calculated ROHF energy in Hartree.
    """
    # Define the molecule: H2+, charge=1, spin=1 (one unpaired electron)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0., 0., -distance/2.)],
        ['H', (0., 0.,  distance/2.)],
    ]
    mol.charge = 1
    mol.spin = 1  # Number of unpaired electrons
    mol.basis = 'cc-pvdz' # Use a standard basis set

    if allow_symmetry_breaking:
        # By setting symmetry to False, we allow the SCF procedure to find a
        # lower-energy solution that may not conform to the molecule's point group.
        mol.symmetry = False
    else:
        # Enforce D2h symmetry (an abelian subgroup of Dinfh). This is the default
        # behavior in many cases and is what leads to the incorrect dissociation.
        mol.symmetry = 'D2h'

    mol.build()

    # Perform a Restricted Open-Shell Hartree-Fock (ROHF) calculation
    mf = scf.ROHF(mol).run(conv_tol=1e-9, quiet=True)

    # For the symmetric case at long distance, we need to ensure we are on the
    # ground state surface (sigma_g) and not the excited state (sigma_u).
    # This can be tricky, but the default guess is usually correct.
    # Forcing a specific orbital occupation can be done if needed.

    return mf.e_tot

def main():
    """
    Main function to calculate and print the potential energy surface.
    """
    # A range of internuclear distances to scan
    distances = np.arange(0.5, 4.1, 0.25)
    # The correct dissociation energy for H + H+ is E(H) = -0.5 Hartree
    correct_dissociation_limit = -0.5

    print("Calculating the H2+ Potential Energy Curve...")
    print("-" * 65)
    print(f"{'Distance (Å)':<15} {'Symmetric E (Ha)':<22} {'Broken-Symmetry E (Ha)':<25}")
    print("-" * 65)

    for r in distances:
        # Calculate energy with symmetry enforced
        energy_symm = calculate_h2_plus_energy(r, allow_symmetry_breaking=False)

        # Calculate energy with symmetry breaking allowed
        energy_broken = calculate_h2_plus_energy(r, allow_symmetry_breaking=True)

        print(f"{r:<15.2f} {energy_symm:<22.6f} {energy_broken:<25.6f}")

    print("-" * 65)
    print(f"Note: The correct dissociation energy limit (R -> inf) is {correct_dissociation_limit:.4f} Hartree.")
    print("\nExplanation:")
    print("The 'Symmetric' calculation incorrectly forces the electron to be delocalized")
    print("over both atoms, leading to a pathologically low energy at long distances.")
    print("\nThe 'Broken-Symmetry' calculation correctly allows the electron to localize on")
    print("one atom, resulting in the correct dissociation behavior where the energy")
    print(f"approaches {correct_dissociation_limit:.4f} Hartree.")


if __name__ == '__main__':
    main()
