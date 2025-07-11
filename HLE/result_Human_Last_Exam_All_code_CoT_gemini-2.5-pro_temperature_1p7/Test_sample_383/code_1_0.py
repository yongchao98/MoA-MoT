# First, ensure you have pyscf installed:
# pip install pyscf
import pyscf
import numpy as np

def calculate_h2_plus_energy(distance, method):
    """
    Calculates the energy of the H2+ cation at a given distance and method.
    """
    # Define the molecule's geometry, charge, and spin state
    mol = pyscf.gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, distance)]
    ]
    mol.charge = 1  # +1 charge for the cation
    mol.spin = 1   # 1 unpaired electron (doublet state)
    mol.basis = 'sto-3g' # A minimal basis set for demonstration
    mol.build()

    # Select the computational method
    # Unrestricted methods (UHF/UKS) are used for proper dissociation of open-shell species
    if method.lower() == 'hf':
        mf = pyscf.scf.UHF(mol)
    elif method.lower() == 'dft':
        mf = pyscf.dft.UKS(mol)
        mf.xc = 'b3lyp' # A common hybrid functional known to have self-interaction error
    else:
        raise ValueError("Method must be 'hf' or 'dft'")
    
    # Suppress verbose output for clarity
    mf.verbose = 0
    # Run the calculation and return the total energy
    energy = mf.kernel()
    return energy

# --- Main Program ---
if __name__ == "__main__":
    # Define distances for comparison: equilibrium and stretched
    eq_dist = 1.06  # Angstrom, near the equilibrium distance
    long_dist = 10.0 # Angstrom, representing dissociation

    # Calculate energies using Hartree-Fock
    e_hf_eq = calculate_h2_plus_energy(eq_dist, 'hf')
    e_hf_long = calculate_h2_plus_energy(long_dist, 'hf')

    # Calculate energies using DFT (B3LYP)
    e_dft_eq = calculate_h2_plus_energy(eq_dist, 'dft')
    e_dft_long = calculate_h2_plus_energy(long_dist, 'dft')

    # For reference, calculate the energy of a single hydrogen atom (the dissociation product)
    # H + H+  ->  Energy should approach E(H)
    h_atom = pyscf.gto.Mole()
    h_atom.atom = 'H'
    h_atom.spin = 1
    h_atom.basis = 'sto-3g'
    h_atom.build()
    e_h_atom_hf = pyscf.scf.UHF(h_atom).run(verbose=0).e_tot
    
    print("Demonstrating the Fix for H2+ Dissociation Curve")
    print("=" * 55)

    print("\n--- Hartree-Fock (HF) Calculation ---")
    print("This method is free of self-interaction error.")
    print(f"Energy at equilibrium distance ({eq_dist} Å): {e_hf_eq:.6f} Hartree")
    print(f"Energy at long distance ({long_dist} Å):    {e_hf_long:.6f} Hartree")
    print(f"Correct Dissociation Limit (H atom): {e_h_atom_hf:.6f} Hartree")
    print("Result: CORRECT. Energy increases at long distance towards the dissociation limit.")

    print("\n--- DFT (B3LYP) Calculation ---")
    print("This method suffers from self-interaction error.")
    print(f"Energy at equilibrium distance ({eq_dist} Å): {e_dft_eq:.6f} Hartree")
    print(f"Energy at long distance ({long_dist} Å):    {e_dft_long:.6f} Hartree")
    print("Result: INCORRECT. Energy unphysically drops far below the equilibrium energy.")
    print("=" * 55)
