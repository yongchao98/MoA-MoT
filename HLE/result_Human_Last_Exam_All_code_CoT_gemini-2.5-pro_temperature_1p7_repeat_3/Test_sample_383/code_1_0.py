import numpy as np
# To run this code, you need to install pyscf: pip install pyscf
from pyscf import gto, scf

def calculate_h2_plus_energy(dist):
    """
    Calculates the energy of H2+ at a given distance using
    both Restricted Open-shell Hartree-Fock (ROHF) and
    a standard DFT functional (LDA).
    """
    # Define the molecule using PySCF's Mole object
    # Charge is +1, and spin is 1 (one unpaired electron, so 2S+1 = 2*0.5 = 1 is incorrect, S=1/2 means spin=1)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, dist)]
    ]
    mol.charge = 1
    mol.spin = 1 # Spin multiplicity is 2S+1. For one electron S=1/2, so 2S+1 = 2. But pyscf uses n_alpha-n_beta, so 1-0=1.
    mol.basis = 'cc-pvdz'
    mol.build()

    # --- Hartree-Fock Calculation (Correct Method) ---
    # ROHF is suitable for open-shell systems. For H2+, it is exact.
    hf_calc = scf.ROHF(mol)
    hf_calc.verbose = 0 # Suppress verbose output
    hf_energy = hf_calc.kernel()

    # --- DFT Calculation (Method with Self-Interaction Error) ---
    # RKS with LDA functional ('lda,vwn') is a common simple DFT.
    dft_calc = scf.RKS(mol)
    dft_calc.xc = 'lda,vwn' # Slater exchange, VWN correlation
    dft_calc.verbose = 0 # Suppress verbose output
    dft_energy = dft_calc.kernel()

    return hf_energy, dft_energy

if __name__ == "__main__":
    # A range of bond distances in Angstroms
    bond_distances = np.arange(0.5, 6.1, 0.25)

    print("Calculating H2+ Potential Energy Curve")
    print("-" * 55)
    print(f"{'Distance (Å)':<15} {'Hartree-Fock (Eh)':<20} {'LDA DFT (Eh)':<20}")
    print("-" * 55)

    for r in bond_distances:
        # Calculate energies
        try:
            e_hf, e_dft = calculate_h2_plus_energy(r)
            # The correct dissociation energy is -0.5 Hartree (Energy of H atom)
            print(f"{r:<15.4f} {e_hf:<20.8f} {e_dft:<20.8f}")
        except Exception as e:
            print(f"{r:<15.4f} {'Calculation failed':<20} {'Calculation failed':<20}")

    print("-" * 55)
    print("\nExplanation:")
    print("The 'Hartree-Fock (Eh)' column shows the correct behavior: energy increases from a minimum and approaches -0.5 Hartree.")
    print("The 'LDA DFT (Eh)' column shows the pathological behavior: the energy incorrectly drops at longer distances due to self-interaction error.")
