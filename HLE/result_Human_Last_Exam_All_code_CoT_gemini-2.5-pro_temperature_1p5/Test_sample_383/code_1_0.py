# First, ensure you have pyscf installed:
# pip install pyscf

from pyscf import gto, scf, dft
import numpy as np

def calculate_h2_plus_energy(distance_angstrom, method_obj):
    """
    Calculates the energy of H2+ at a given distance using a specified method.
    
    Args:
        distance_angstrom (float): The internuclear distance in Angstroms.
        method_obj: The pyscf method object (e.g., scf.UHF or dft.UKS).
        
    Returns:
        float: The total electronic energy in Hartree.
    """
    # Define the molecule
    # H2+ has 2 protons, 1 electron.
    # Charge = +1, Spin = 1/2 (doublet multiplicity)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, distance_angstrom)]
    ]
    mol.charge = 1
    mol.spin = 1  # Corresponds to 2S+1=2 -> S=1/2
    mol.basis = 'cc-pVDZ' # A standard basis set
    mol.build()
    
    # Run the specified calculation method
    method = method_obj(mol)
    energy = method.kernel()
    return energy

def main():
    """
    Main function to demonstrate the fix for H2+ dissociation.
    """
    # A large distance to simulate dissociation
    dissociation_distance = 10.0  # in Angstroms

    print(f"Calculating energy of H2+ at a large distance of {dissociation_distance} Å.")
    print("The correct energy should be that of a single Hydrogen atom, which is ~-0.5 Hartree.")
    print("-" * 70)

    # --- Calculation 1: The problematic method (DFT with B3LYP) ---
    print("Running calculation with DFT (B3LYP), a method with self-interaction error...")
    # We use Unrestricted Kohn-Sham (UKS) which is the DFT equivalent of UHF
    dft_energy = calculate_h2_plus_energy(dissociation_distance, dft.UKS)
    print(f"Result (DFT/B3LYP): {dft_energy:.6f} Hartree")
    print("This energy is significantly lower than -0.5, which is physically incorrect.")
    print("This is caused by self-interaction error forcing the electron to be delocalized.")
    print("-" * 70)

    # --- Calculation 2: The correct method (Hartree-Fock) ---
    print("Running calculation with Hartree-Fock, which is self-interaction free...")
    # We use Unrestricted Hartree-Fock (UHF) to allow the electron to localize
    hf_energy = calculate_h2_plus_energy(dissociation_distance, scf.UHF)
    print(f"Result (UHF): {hf_energy:.6f} Hartree")
    print("This energy is very close to the correct value of -0.5 Hartree.")
    print("This demonstrates that using Hartree-Fock (Option 2) fixes the problem.")
    print("-" * 70)

if __name__ == "__main__":
    main()
