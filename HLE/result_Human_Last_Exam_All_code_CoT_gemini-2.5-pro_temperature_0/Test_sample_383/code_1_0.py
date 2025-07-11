import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_pes():
    """
    Calculates and prints the potential energy surface for H2+ using
    both a DFT functional (PBE) and Hartree-Fock (UHF) to demonstrate
    the correction of the self-interaction error.
    """
    # --- Explanation ---
    # The problem described is a classic example of self-interaction error (SIE)
    # in many common density functional approximations (DFAs). This error causes
    # an unphysical lowering of the energy as the H2+ bond is stretched.
    #
    # Statement 3 correctly identifies the type of system where this occurs.
    # Statement 2 provides the solution: Hartree-Fock (HF) theory is self-interaction
    # free for one-electron systems. Using Unrestricted Hartree-Fock (UHF)
    # correctly models the dissociation.

    # Define a range of internuclear distances in Angstroms
    distances = np.arange(0.6, 5.1, 0.2)

    print("Calculating Potential Energy Curves for H2+")
    print("-" * 60)
    # The final "equation" is Energy(Distance). We print each component.
    print(f"{'Distance (A)':<15} {'PBE DFT Energy (Ha)':<22} {'UHF Energy (Ha)':<20}")
    print("-" * 60)

    # Loop over each distance to calculate the energy
    for r in distances:
        # Define the H2+ molecule in PySCF
        # charge=1, spin=1 (one unpaired electron)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, r)]
        ]
        mol.charge = 1
        mol.spin = 1
        mol.basis = 'cc-pvdz'
        # Suppress PySCF's default output for cleaner printing
        mol.verbose = 0
        mol.build()

        # --- Problematic Calculation: DFT with PBE functional ---
        # We use Unrestricted Kohn-Sham (UKS), the open-shell version of DFT
        mf_dft = dft.UKS(mol)
        mf_dft.xc = 'pbe,pbe' # PBE functional
        e_dft = mf_dft.kernel()

        # --- Fixed Calculation: Unrestricted Hartree-Fock (UHF) ---
        mf_hf = scf.UHF(mol)
        e_hf = mf_hf.kernel()

        # Output each number (distance, energy1, energy2) for the final result
        print(f"{r:<15.2f} {e_dft:<22.6f} {e_hf:<20.6f}")

    print("-" * 60)
    print("\nAnalysis:")
    print("Note how the PBE DFT energy incorrectly drops at long distances.")
    print("The UHF energy correctly flattens out near -0.5 Hartree, the")
    print("correct energy for a separated H atom and a proton.")

if __name__ == '__main__':
    calculate_h2_plus_pes()