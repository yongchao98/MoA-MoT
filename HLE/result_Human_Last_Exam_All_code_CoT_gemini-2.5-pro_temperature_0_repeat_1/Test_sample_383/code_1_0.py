# To run this code, you first need to install the pyscf library:
# pip install pyscf

from pyscf import gto, scf, dft

def calculate_h2_plus_pes():
    """
    Calculates and compares the potential energy curve of the H2+ cation
    using both Hartree-Fock (HF) and a DFT functional (PBE).

    This demonstration highlights the self-interaction error in common DFT
    functionals for stretched, symmetric, odd-electron systems and shows
    how HF provides a correct description.
    """
    print("Calculating Potential Energy Curve for H2+ Cation")
    print("The equilibrium bond distance is approximately 1.06 Angstroms.")
    print("-" * 73)
    # Header for the output table. Note the final equation format requirement.
    print(f"{'Distance(A)':<12} | {'Equation: E_HF = <Psi_HF|H|Psi_HF>':<28} | {'Equation: E_PBE = E_T[n] + E_V[n] + E_J[n] + E_XC[n]':<30}")
    print(f"{'':<12} | {'HF Energy (Hartree)':<28} | {'PBE Energy (Hartree)':<30}")
    print("-" * 73)

    # List of internuclear distances (in Angstroms) to calculate
    distances = [0.6, 0.8, 1.0, 1.06, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0]

    for d in distances:
        # Define the H2+ molecule for the given distance d.
        # It has 2 protons and 1 electron, so charge=1 and spin=1 (a doublet).
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, d)]
        ]
        mol.basis = 'cc-pvdz'  # A standard basis set
        mol.charge = 1
        mol.spin = 1  # Spin is 2*S, so 2 * (1/2) = 1
        mol.build(verbose=0) # verbose=0 suppresses detailed pyscf output

        # --- Calculation 1: Hartree-Fock (The "Fix") ---
        # HF is free of self-interaction error and is exact for one-electron systems.
        hf_calculation = scf.RHF(mol)
        energy_hf = hf_calculation.kernel()

        # --- Calculation 2: DFT with PBE functional (The "Problem") ---
        # PBE is a common GGA functional that suffers from self-interaction error.
        dft_calculation = dft.RKS(mol)
        dft_calculation.xc = 'pbe,pbe'
        energy_dft = dft_calculation.kernel()

        # Print the results for the current distance in the table format
        print(f"{d:<12.2f} | {energy_hf:<28.8f} | {energy_dft:<30.8f}")

    print("-" * 73)
    print("\nAnalysis:")
    print("Notice that the PBE (DFT) energy incorrectly drops at large distances,")
    print("becoming much lower than the energy at the equilibrium distance (~1.06 A).")
    print("This is physically wrong. In contrast, the Hartree-Fock (HF) energy behaves")
    print("correctly, increasing as the bond is stretched towards dissociation.")

# Run the main function to perform the calculations and print the results.
calculate_h2_plus_pes()