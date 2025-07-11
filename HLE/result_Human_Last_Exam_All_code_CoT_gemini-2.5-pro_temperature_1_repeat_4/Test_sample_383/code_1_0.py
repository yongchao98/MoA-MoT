import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_pes():
    """
    Calculates the potential energy surface for H2+ using
    Restricted Open-Shell DFT (which fails) and
    Restricted Open-Shell Hartree-Fock (which succeeds).

    This demonstrates the self-interaction error problem in some DFT functionals
    for stretched radical cations.
    """
    # Define a range of internuclear distances in Angstroms
    distances = np.arange(0.6, 5.1, 0.2)

    print("H2+ Potential Energy Curve Calculation")
    print("-" * 65)
    print("{:<15} {:<25} {:<25}".format("Distance (A)", "ROKS/B3LYP Energy (Ha)", "ROHF Energy (Ha)"))
    print("-" * 65)

    # The exact energy of a hydrogen atom (the dissociation limit) is -0.5 Ha
    dissociation_limit = -0.5

    for r in distances:
        # Define the molecule at the given distance
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, r)]]
        mol.charge = 1  # H2+ cation
        mol.spin = 1    # One unpaired electron (Doublet state)
        mol.basis = 'cc-pvdz' # Use a decent basis set
        mol.build()

        # --- Failing Method: Restricted Open-Shell DFT (ROKS) with B3LYP ---
        # ROKS enforces symmetry, and B3LYP suffers from self-interaction error,
        # leading to an unphysical energy drop at large distances.
        # We suppress verbose output with mol.verbose=0
        mol.verbose = 0
        mf_roks = scf.ROKS(mol)
        mf_roks.xc = 'b3lyp'
        e_roks = mf_roks.kernel()

        # --- Correct Method: Restricted Open-Shell Hartree-Fock (ROHF) ---
        # ROHF is free of self-interaction error for a one-electron system
        # and gives the correct dissociation behavior.
        mf_rohf = scf.ROHF(mol)
        e_rohf = mf_rohf.kernel()

        # Print the results for the current distance
        print("{:<15.2f} {:<25.8f} {:<25.8f}".format(r, e_roks, e_rohf))

    print("-" * 65)
    print(f"Note: The correct dissociation limit is E(H) + E(H+) = {dissociation_limit:.1f} Hartree.")
    print("Observe how the ROKS/B3LYP energy incorrectly drops far below the")
    print("minimum, while the ROHF energy correctly approaches the -0.5 Ha limit.")


if __name__ == '__main__':
    # To run this code, you need to install pyscf:
    # pip install pyscf
    calculate_h2_plus_pes()

<<<G>>>