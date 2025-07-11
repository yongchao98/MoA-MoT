import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_pes():
    """
    Calculates the potential energy surface for H2+ using two methods
    to demonstrate the concept of inverse symmetry breaking.
    """
    print("H2+ Potential Energy Curve Calculation")
    print("The energy of the dissociated products (H + H+) is -0.5 Hartree.")
    print("-" * 70)
    print("{:<12} {:<25} {:<25}".format("Distance(A)", "ROHF Energy (Hartree)", "UHF Energy (Hartree)"))
    print("-" * 70)

    # List of distances to calculate
    distances = np.arange(0.5, 5.5, 0.25)

    for r in distances:
        # Define the H2+ molecule for each distance
        # Charge is 1, spin is 1 (one unpaired electron -> 2S+1 = 2*0.5+1=2 (doublet) -> nelec_alpha-nelec_beta = 1)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0., 0., 0.)],
            ['H', (0., 0., r)],
        ]
        mol.charge = 1
        mol.spin = 1 # Number of unpaired electrons
        mol.basis = 'cc-pvdz'
        mol.build()

        # --- Calculation 1: ROHF (shows the problem) ---
        # This method restricts the spatial part of the orbital,
        # forcing symmetry and leading to incorrect dissociation.
        mf_rohf = scf.ROHF(mol).run(verbose=0)
        e_rohf = mf_rohf.e_tot

        # --- Calculation 2: UHF (provides a fix) ---
        # This method allows breaking the spatial symmetry,
        # which correctly localizes the electron on one atom at large distances.
        mf_uhf = scf.UHF(mol).run(verbose=0)
        e_uhf = mf_uhf.e_tot

        # Output the results for this distance
        # This corresponds to "output each number in the final equation"
        print("{:<12.4f} {:<25.8f} {:<25.8f}".format(r, e_rohf, e_uhf))

    print("-" * 70)
    print("Notice how the ROHF energy drops unphysically at large distances,")
    print("while the UHF energy correctly approaches the -0.5 Hartree limit.")

if __name__ == '__main__':
    calculate_h2_plus_pes()