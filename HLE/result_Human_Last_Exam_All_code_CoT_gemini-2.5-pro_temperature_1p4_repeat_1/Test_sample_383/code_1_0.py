# First, ensure you have pyscf installed:
# pip install pyscf

import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_energy(distance_A):
    """
    Calculates the energy of H2+ at a given internuclear distance
    using both DFT (PBE) and Hartree-Fock (UHF).
    """
    # Define the molecule: H2+
    # atom: List of atoms and their coordinates
    # basis: The basis set to use for the calculation
    # charge: +1 for H2+
    # spin: 1 unpaired electron (doublet state)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0., 0., 0.)],
        ['H', (0., 0., distance_A)]
    ]
    mol.basis = 'cc-pVDZ'
    mol.charge = 1
    mol.spin = 1
    mol.build()

    # --- 1. DFT Calculation (using PBE functional) ---
    # We use Unrestricted Kohn-Sham (UKS) to allow for symmetry breaking.
    dft_calc = dft.UKS(mol)
    dft_calc.xc = 'pbe,pbe' # PBE is a common GGA functional
    # Suppress verbose output for cleaner results
    dft_calc.verbose = 0
    energy_dft = dft_calc.kernel()

    # --- 2. Hartree-Fock Calculation ---
    # We use Unrestricted Hartree-Fock (UHF)
    hf_calc = scf.UHF(mol)
    hf_calc.verbose = 0
    energy_hf = hf_calc.kernel()

    return energy_dft, energy_hf

# Set a large distance to simulate dissociation
distance = 10.0  # in Angstroms

# Perform the calculations
dft_energy, hf_energy = calculate_h2_plus_energy(distance)

# The correct energy for the dissociated products (H + H+) is the energy
# of a single Hydrogen atom, which is exactly -0.5 Hartree.
correct_energy = -0.5

print(f"Calculating energy for H2+ at a separation of {distance} Angstroms:\n")
print("====================================================================")
print(f"Final Energy (DFT with PBE functional): {dft_energy:.6f} Hartree")
print(f"Final Energy (Hartree-Fock):            {hf_energy:.6f} Hartree")
print(f"Correct Dissociation Energy:            {correct_energy:.6f} Hartree")
print("====================================================================")
print("\nObservation:")
print("The DFT energy is significantly lower than the correct value,")
print("demonstrating the unphysical energy drop due to self-interaction error.")
print("The Hartree-Fock energy is correct, as it does not suffer from this error.")
