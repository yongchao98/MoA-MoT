import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_pes():
    """
    Calculates the potential energy surface (PES) for the H2+ molecular cation
    using both Density Functional Theory (PBE) and Hartree-Fock (HF).

    This demonstrates the self-interaction error (SIE) in the PBE functional,
    which leads to an incorrect dissociation curve, and shows how HF, being
    exact for a one-electron system, provides the correct curve.
    """

    # Define a range of internuclear distances in Angstroms
    distances = np.arange(0.6, 5.1, 0.2)
    
    # Lists to store the calculated energies
    hf_energies = []
    dft_energies = []

    print("Calculating Potential Energy Curve for H2+...")
    print("-" * 55)
    print(f"{'Distance (A)':<15} | {'PBE Energy (Hartree)':<20} | {'HF Energy (Hartree)':<20}")
    print("-" * 55)

    for r in distances:
        # Define the H2+ molecule for the current distance
        # pyscf uses Angstroms for geometry definition
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, r)]
        ]
        mol.basis = 'cc-pvdz'  # A decent quality basis set
        mol.charge = 1           # H2+ has a +1 charge
        mol.spin = 1             # One unpaired electron (doublet state)
        mol.build()

        # 1. DFT calculation with PBE functional (shows the error)
        # Using RKS (Restricted Kohn-Sham) for this doublet system
        dft_calc = dft.RKS(mol)
        dft_calc.xc = 'pbe'  # PBE is a common GGA functional with known SIE
        dft_calc.verbose = 0 # Suppress detailed output for cleaner printing
        e_dft = dft_calc.kernel()
        dft_energies.append(e_dft)

        # 2. Hartree-Fock calculation (the correct method for this system)
        # Using RHF (Restricted Hartree-Fock)
        hf_calc = scf.RHF(mol)
        hf_calc.verbose = 0 # Suppress detailed output
        e_hf = hf_calc.kernel()
        hf_energies.append(e_hf)
        
        # Print the results for each distance as they are calculated
        print(f"{r:<15.2f} | {e_dft:<20.6f} | {e_hf:<20.6f}")

    print("-" * 55)
    print("\nCalculation complete. Observe the results:")
    print("1. The HF energy starts at a minimum and correctly increases towards a plateau.")
    print("2. The PBE (DFT) energy incorrectly drops at longer distances due to Self-Interaction Error.")

# Execute the function
calculate_h2_plus_pes()
