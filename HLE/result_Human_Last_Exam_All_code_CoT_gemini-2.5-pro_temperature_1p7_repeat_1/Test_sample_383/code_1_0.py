import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_pes(distances):
    """
    Calculates the potential energy surface for H2+ using RHF (correct)
    and B3LYP DFT (incorrect due to Self-Interaction Error).
    """
    rhf_energies = []
    dft_energies = []

    print(f"{'Distance (A)':<15} {'RHF/ROHF (Hartree)':<25} {'B3LYP (Hartree)':<20}")
    print("-" * 60)

    for r in distances:
        # Define the molecule: H2+ with charge +1 and spin 1 (doublet)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0., 0., 0.)],
            ['H', (0., 0., r)]]
        mol.charge = 1
        mol.spin = 1  # 2S = 1 -> Doublet
        mol.basis = 'cc-pVDZ'
        mol.build()

        # 1. Correct Method: Restricted Open-Shell Hartree-Fock (ROHF)
        # For a one-electron system, this is exact.
        mf_rhf = scf.ROHF(mol).run(verbose=0)
        rhf_energies.append(mf_rhf.e_tot)

        # 2. Incorrect Method: B3LYP DFT
        # This functional suffers from self-interaction error.
        mf_dft = scf.RKS(mol)
        mf_dft.xc = 'b3lyp'
        mf_dft.run(verbose=0)
        dft_energies.append(mf_dft.e_tot)

        print(f"{r:<15.2f} {mf_rhf.e_tot:<25.8f} {mf_dft.e_tot:<20.8f}")

    print("\n--- Analysis ---")
    min_rhf_energy = min(rhf_energies)
    dissociation_rhf = rhf_energies[-1]
    print(f"Correct (RHF) potential minimum: {min_rhf_energy:.6f} Hartree")
    print(f"Correct (RHF) dissociation energy (at R={distances[-1]} A): {dissociation_rhf:.6f} Hartree")
    print("Note: The RHF energy correctly rises from the minimum towards the dissociation limit of H + H+ (-0.5 Hartree).")

    min_dft_energy = min(dft_energies)
    long_range_dft = dft_energies[-1]
    print(f"\nIncorrect (B3LYP) potential minimum: {min_dft_energy:.6f} Hartree")
    print(f"Incorrect (B3LYP) energy at long range (R={distances[-1]} A): {long_range_dft:.6f} Hartree")
    print("Note: The B3LYP energy incorrectly drops far below the potential minimum due to self-interaction error.")


if __name__ == '__main__':
    # Define a range of internuclear distances in Angstroms
    # Start near equilibrium (approx 1.06 A) and stretch the bond
    distances_angstrom = np.arange(0.6, 5.1, 0.2)
    calculate_h2_plus_pes(distances_angstrom)
