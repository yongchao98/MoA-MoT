# The user needs to have psi4 installed. It can be installed via conda or pip:
# conda install psi4 -c psi4
# or
# pip install psi4

import psi4
import numpy as np

def calculate_h2_plus_energy():
    """
    Demonstrates the correct and incorrect calculation of the H2+ potential energy curve.
    """
    print("This script demonstrates the 'inverse symmetry breaking' or 'self-interaction error' problem for H2+ and its solution.")
    print("We calculate energies at various distances using a problematic method (UKS/B3LYP) and the correct method (UHF).\n")
    
    # Define the molecular structure template and distances
    mol_template = """
    1 2
    H 0 0 0
    H 0 0 {R}
    symmetry c1
    """
    distances = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]

    # --- Calculation with a problematic method (UKS/B3LYP) ---
    # This method suffers from self-interaction error, causing an unphysical energy drop.
    psi4.core.clean()
    psi4.set_options({
        'basis': 'cc-pvdz',
        'scf_type': 'uks'  # Unrestricted Kohn-Sham for open-shell DFT
    })
    
    b3lyp_energies = []
    for r in distances:
        mol = psi4.geometry(mol_template.format(R=r))
        # quiet=True suppresses detailed psi4 output for each step
        energy = psi4.energy('b3lyp', molecule=mol, quiet=True)
        b3lyp_energies.append(energy)
        
    # --- Calculation with the correct method (UHF) ---
    # Unrestricted Hartree-Fock (UHF) is self-interaction free and allows for correct
    # symmetry breaking as the bond dissociates.
    psi4.core.clean()
    psi4.set_options({
        'basis': 'cc-pvdz',
        'scf_type': 'uhf'  # Unrestricted Hartree-Fock
    })
    
    uhf_energies = []
    for r in distances:
        mol = psi4.geometry(mol_template.format(R=r))
        energy = psi4.energy('hf', molecule=mol, quiet=True)
        uhf_energies.append(energy)

    # --- Print the results ---
    print("--- Comparison of Calculated Energies (in Hartree) ---")
    print(f"{'Distance (A)':<15} {'UKS/B3LYP (Error)':<20} {'UHF (Correct)':<20}")
    print("-" * 60)
    for i, r in enumerate(distances):
        print(f"{r:<15.1f} {b3lyp_energies[i]:<20.6f} {uhf_energies[i]:<20.6f}")
        
    print("\n--- Analysis ---")
    print("Observe the B3LYP energy at R >= 3.0 A. It drops unphysically low, below the equilibrium energy (at R=1.0 A).")
    print("In contrast, the UHF energy correctly rises from equilibrium and levels off towards the correct dissociation limit of -0.5 Hartree.")

calculate_h2_plus_energy()