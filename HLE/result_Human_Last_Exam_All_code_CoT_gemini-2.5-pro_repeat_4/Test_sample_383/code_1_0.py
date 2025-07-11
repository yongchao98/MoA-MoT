import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_energy(distance, method):
    """
    Calculates the energy of H2+ at a given distance using a specified method.
    
    Args:
        distance (float): The internuclear distance in Angstroms.
        method (str): 'hf' for Hartree-Fock or 'dft' for DFT (LDA).
        
    Returns:
        float: The total energy in Hartrees.
    """
    # Define the molecule: H2+
    # Charge = +1, Spin = 1 (doublet, 1 unpaired electron)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, distance)]
    ]
    mol.charge = 1
    mol.spin = 1
    mol.basis = 'cc-pVDZ' # Use a standard basis set
    mol.build()
    
    # Select the calculation method
    if method.lower() == 'hf':
        # Use Restricted Open-shell Hartree-Fock (ROHF)
        calc = scf.ROHF(mol)
    elif method.lower() == 'dft':
        # Use Restricted Open-shell Kohn-Sham DFT with LDA functional
        calc = scf.ROKS(mol)
        calc.xc = 'lda,vwn' # A simple DFT functional known to have large SIE
    else:
        raise ValueError("Method must be 'hf' or 'dft'")

    # Run the calculation and suppress verbose output
    calc.verbose = 0
    energy = calc.kernel()
    return energy

def main():
    """
    Main function to calculate and print the potential energy curves.
    """
    # A range of bond distances to scan, from near equilibrium to dissociated
    distances = np.arange(0.5, 6.1, 0.25)
    
    print("Calculating Potential Energy Curve for H2+")
    print("="*50)
    print(f"{'Distance (A)':<15} {'E(HF) / Hartree':<20} {'E(DFT/LDA) / Hartree':<20}")
    print("-" * 50)
    
    for d in distances:
        # Calculate energy using Hartree-Fock (the correct method here)
        energy_hf = calculate_h2_plus_energy(d, 'hf')
        
        # Calculate energy using DFT (the problematic method)
        energy_dft = calculate_h2_plus_energy(d, 'dft')
        
        # Print each number in the results clearly
        print(f"{d:<15.2f} {energy_hf:<20.8f} {energy_dft:<20.8f}")
        
    print("-" * 50)
    print("\nObservation:")
    print("The E(HF) column shows the correct behavior: energy increases from the minimum")
    print("and flattens out around -0.5 Hartree (the energy of a Hydrogen atom).")
    print("The E(DFT/LDA) column shows the incorrect behavior: the energy keeps dropping")
    print("unphysically at larger distances due to self-interaction error.")

if __name__ == "__main__":
    main()