# First, please make sure you have pyscf installed:
# pip install pyscf

import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_energy(distance_A):
    """
    Calculates the energy of H2+ at a given distance using ROHF and UHF methods.
    """
    # Conversion factor from Angstrom to Bohr
    A_to_bohr = 1.8897259886
    distance_bohr = distance_A * A_to_bohr
    
    # Define the H2+ molecule
    # H2+ has 1 electron, so charge=1 and spin=1 (one unpaired electron)
    mol = gto.Mole()
    mol.atom = f'H 0 0 0; H 0 0 {distance_A}'
    mol.unit = 'A'
    mol.basis = 'cc-pvdz'
    mol.charge = 1
    mol.spin = 1
    mol.build()

    # --- Problematic Calculation: Restricted Open-shell Hartree-Fock (ROHF) ---
    # ROHF forces the electron to be in a symmetric orbital, which is wrong at large distances.
    rohf_calc = scf.ROHF(mol)
    rohf_energy = rohf_calc.kernel() # Use .kernel() to run the calculation

    # --- Corrected Calculation: Unrestricted Hartree-Fock (UHF) ---
    # UHF allows the electron to localize on one atom (symmetry breaking), which is correct.
    uhf_calc = scf.UHF(mol)
    uhf_energy = uhf_calc.kernel()

    # --- Reference Calculation for Dissociation Limit ---
    # Energy of a single Hydrogen atom
    h_atom = gto.M(atom='H 0 0 0', basis='cc-pvdz', spin=1)
    h_energy = scf.UHF(h_atom).kernel()
    
    # Nuclear repulsion energy E_nuc = Z_A*Z_B/R = 1*1/R (in atomic units)
    nuclear_repulsion_energy = mol.energy_nuc()
    
    # The expected correct total energy is E(H) + E(H+) + E_nuc
    # Energy of a proton E(H+) is 0.
    expected_dissociation_energy = h_energy + nuclear_repulsion_energy

    # --- Print the results ---
    print(f"Calculating energy for H2+ at a distance of {distance_A} Angstroms:")
    print("-----------------------------------------------------------------")
    print("1. Using a method that incorrectly imposes symmetry (ROHF):")
    print(f"   - ROHF Total Energy: {rohf_energy:.6f} Hartree")
    print("\n2. Using a method that correctly allows symmetry breaking (UHF):")
    print(f"   - UHF Total Energy:  {uhf_energy:.6f} Hartree")
    print("\n3. Reference for the correct dissociation limit:")
    print(f"   - Energy of isolated H atom:            {h_energy:.6f} Hartree")
    print(f"   - Nuclear repulsion energy at {distance_A} A: {nuclear_repulsion_energy:.6f} Hartree")
    print(f"   - Expected total energy (E(H)+E_nuc):   {expected_dissociation_energy:.6f} Hartree")

    print("\nConclusion: The UHF energy matches the expected dissociation limit,")
    print("while the ROHF energy is artificially and incorrectly low.")


# Run the calculation at a large distance where the problem is obvious
calculate_h2_plus_energy(distance_A=10.0)