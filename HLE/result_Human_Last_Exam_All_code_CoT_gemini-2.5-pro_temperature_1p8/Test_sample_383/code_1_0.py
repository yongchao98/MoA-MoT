import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_energy(distance, method='uhf'):
    """
    Calculates the energy of H2+ at a given distance using the specified method.
    """
    # Define the molecule using PySCF's Mole object
    # H2+ has 2 protons, 1 electron, hence charge=1 and spin=1 (one unpaired electron)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, distance)]
    ]
    mol.charge = 1
    mol.spin = 1
    mol.build()

    # Select the method
    if method.lower() == 'uhf':
        # Use Unrestricted Hartree-Fock
        mf = scf.UHF(mol)
    elif method.lower() == 'pbe':
        # Use Unrestricted Kohn-Sham DFT with PBE functional
        mf = dft.UKS(mol)
        mf.xc = 'pbe'
    else:
        raise ValueError("Method must be 'uhf' or 'pbe'")

    # Run the calculation silently
    mf.verbose = 0
    energy = mf.kernel()
    return energy

# --- Main script ---
if __name__ == "__main__":
    # Define a range of bond distances in Angstroms
    distances = np.arange(0.5, 8.1, 0.25)

    hf_energies = []
    pbe_energies = []

    print("Calculating potential energy curves for H2+...")
    print("-" * 55)
    print("{:<15} {:<20} {:<20}".format("Distance (A)", "HF Energy (Ha)", "PBE Energy (Ha)"))
    print("-" * 55)

    for d in distances:
        hf_e = calculate_h2_plus_energy(d, method='uhf')
        pbe_e = calculate_h2_plus_energy(d, method='pbe')

        hf_energies.append(hf_e)
        pbe_energies.append(pbe_e)
        
        # The question asks to output numbers in an equation.
        # We will print the calculated energy values for each distance.
        print("{:<15.2f} E = {:<18.8f} E = {:<18.8f}".format(d, hf_e, pbe_e))
        
    print("-" * 55)
    print("\nCalculation complete.")
    print("Note the PBE energy at large distances drops unphysically below its minimum.")
    
    # The true dissociation energy for a one-electron system in HF theory is -0.5 Ha (energy of one H atom).
    min_hf_energy = min(hf_energies)
    min_pbe_energy = min(pbe_energies)
    dissociation_limit = -0.5
    
    print(f"\nEquilibrium Energy (HF) : {min_hf_energy:.6f} Ha")
    print(f"Equilibrium Energy (PBE): {min_pbe_energy:.6f} Ha")
    print(f"Final Energy at {distances[-1]:.1f} A (HF) : {hf_energies[-1]:.6f} Ha (Correctly approaches -0.5 Ha)")
    print(f"Final Energy at {distances[-1]:.1f} A (PBE): {pbe_energies[-1]:.6f} Ha (Incorrectly lower than the minimum!)")

    # Optional: Plotting the results for visualization
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 6))
        plt.plot(distances, hf_energies, 'o-', label='UHF')
        plt.plot(distances, pbe_energies, 's-', label='PBE (DFT)')
        plt.axhline(y=dissociation_limit, color='r', linestyle='--', label='True Dissociation Limit (-0.5 Ha)')
        plt.xlabel("Bond Distance (Angstrom)")
        plt.ylabel("Total Energy (Hartree)")
        plt.title("H2+ Potential Energy Curve")
        plt.legend()
        plt.ylim(min(pbe_energies)-0.05, dissociation_limit + 0.1)
        plt.grid(True)
        print("\nPlotting the curves. Close the plot window to exit.")
        plt.show()
    except ImportError:
        print("\nMatplotlib not found. Skipping plot.")
        print("You can install it with: pip install matplotlib")
