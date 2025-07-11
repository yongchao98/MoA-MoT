import psi4
import numpy as np

# --- Settings ---
# Set memory and number of threads
psi4.set_memory('500 MB')
psi4.set_num_threads(2)

# Set the basis set
psi4.set_options({'basis': 'cc-pvdz'})

# --- H2+ Molecule Definition ---
# H2+ has one electron, so it's a doublet (multiplicity 2) and has a +1 charge.
h2_cation = """
{dist} 1
H
H 1 R
symmetry c1
charge = 1
multiplicity = 2
"""

# --- Calculation ---
# Define the range of distances to scan (in Angstroms)
distances = np.arange(0.7, 8.1, 0.2)

print("Calculating Potential Energy Surface for H2+")
print("Dissociation Limit (Energy of H atom): -0.5 Hartree\n")
print(f"{'Distance (A)':<15} {'HF Energy (Ha)':<20} {'PBE Energy (Ha)':<20}")
print("-" * 55)

for dist in distances:
    # Update the geometry for the current distance
    mol = psi4.geometry(h2_cation.format(dist=dist))

    # Calculate energy using Hartree-Fock (correct method for this case)
    try:
        energy_hf = psi4.energy('hf', molecule=mol)
    except psi4.SCFConvergenceError as e:
        energy_hf = float('nan')
        print(f"HF failed to converge at R={dist:.2f} A")


    # Calculate energy using PBE DFT functional (shows self-interaction error)
    try:
        energy_pbe = psi4.energy('pbe', molecule=mol)
    except psi4.SCFConvergenceError as e:
        energy_pbe = float('nan')
        print(f"PBE failed to converge at R={dist:.2f} A")


    # Print the results for the current distance
    print(f"{dist:<15.2f} {energy_hf:<20.8f} {energy_pbe:<20.8f}")
