import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_pec():
    """
    Calculates the Potential Energy Curve (PEC) for H2+
    using both Hartree-Fock (HF) and DFT (PBE functional).

    This demonstrates the failure of standard DFT for this system,
    where the energy spuriously drops at large internuclear distances,
    and how HF provides the correct behavior.
    """
    # Define a range of internuclear distances to scan
    distances = np.arange(0.5, 6.1, 0.25)

    print("Calculating H2+ Potential Energy Curve")
    print("-" * 50)
    print("{:<10s} {:<20s} {:<20s}".format("R (Ang)", "E(HF) [Hartree]", "E(DFT-PBE) [Hartree]"))
    print("-" * 50)

    hf_energies = []
    dft_energies = []

    for r in distances:
        # Define the H2+ molecule for the current distance r
        # We set charge=1 and spin=1 (doublet state for one electron)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, -r / 2)],
            ['H', (0, 0, r / 2)],
        ]
        mol.basis = 'cc-pVDZ'
        mol.charge = 1
        mol.spin = 1
        mol.build()

        # 1. Hartree-Fock Calculation (the correct method here)
        mf_hf = scf.RHF(mol).run(conv_tol=1e-9, verbose=0)
        hf_energies.append(mf_hf.e_tot)

        # 2. DFT Calculation with a standard GGA functional (PBE)
        # This method will show the pathological behavior.
        mf_dft = dft.RKS(mol)
        mf_dft.xc = 'pbe,pbe' # PBE functional
        mf_dft.conv_tol = 1e-9
        # Suppress convergence warnings that may occur at large distances
        mf_dft.verbose = 0
        mf_dft.run()
        dft_energies.append(mf_dft.e_tot)

        # Print the energy at each point in the final equation format
        print("At R = {:.2f} Angstrom:".format(r))
        print("  Equation (HF): E_total = {:.8f} Hartree".format(mf_hf.e_tot))
        print("  Equation (DFT): E_total = {:.8f} Hartree".format(mf_dft.e_tot))
        print()


    print("\n--- Summary ---")
    print("{:<10s} {:<20s} {:<20s}".format("R (Ang)", "E(HF) [Hartree]", "E(DFT-PBE) [Hartree]"))
    print("-" * 50)
    for r, ehf, edft in zip(distances, hf_energies, dft_energies):
        print("{:<10.2f} {:<20.8f} {:<20.8f}".format(r, ehf, edft))

    print("\nNote: Observe how the DFT energy drops unphysically at large R,")
    print("while the HF energy correctly rises and flattens out.")

if __name__ == '__main__':
    # To run this code, you need to have pyscf and numpy installed:
    # pip install pyscf numpy
    calculate_h2_plus_pec()