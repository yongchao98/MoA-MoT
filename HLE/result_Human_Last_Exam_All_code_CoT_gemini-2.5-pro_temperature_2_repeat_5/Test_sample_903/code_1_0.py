#
# This script determines the coordinated atoms in a chemical reaction based on established coordination chemistry principles.
#

# Step 1: Define the reactants and their potential coordinating atoms.

# The ligand is 1,3-Di[3-(2’-pyridyl)pyrazol-1-ylmethyl]-2,4,6-triethylbenzene.
# Each of the two arms has one pyridyl nitrogen and one pyrazolyl nitrogen that can coordinate.
ligand_donors = {'N': 4}

# The metal salt is ZnBr2.
metal_center = 'Zn'
anion_donors = {'Br': 2}

# The solvent is methanol, a potential but weak oxygen donor.
solvent_donors = {'O': 1} # per molecule

print("Analysis of Reactants:")
print(f"Ligand provides {ligand_donors['N']} Nitrogen donor atoms.")
print(f"Metal salt provides the '{metal_center}' center and {anion_donors['Br']} Bromide ions.")
print("---")


# Step 2: Apply coordination principles.
# The strong, tetradentate ligand will coordinate first due to the chelate effect.
coordinated_atoms = []
for atom, count in ligand_donors.items():
    coordinated_atoms.extend([atom] * count)

print("Step-by-step Coordination:")
print(f"First, the four nitrogen atoms from the tetradentate ligand coordinate to the Zinc center.")
print(f"Coordinated atoms so far: {', '.join(coordinated_atoms)}")
current_coordination_number = len(coordinated_atoms)
print(f"Current coordination number: {current_coordination_number}")
print("---")


# Step 3: Check for further coordination to achieve a stable geometry.
# Zn(II) commonly forms 6-coordinate (octahedral) complexes.
# The bromide ions are good ligands and are available to fill the remaining sites.
print("Zinc(II) commonly forms 6-coordinate complexes.")
print("The two available Bromide ions are good ligands and will coordinate.")

for atom, count in anion_donors.items():
    coordinated_atoms.extend([atom] * count)
    
final_coordination_number = len(coordinated_atoms)

# Sort the list for consistent representation
coordinated_atoms.sort()

# Final output
print("Final coordinated atoms after Bromides bind:")
# The final result prints each atom to represent the final coordinated structure.
print(f"The final coordinated atoms are: {coordinated_atoms[0]}, {coordinated_atoms[1]}, {coordinated_atoms[2]}, {coordinated_atoms[3]}, {coordinated_atoms[4]}, {coordinated_atoms[5]}")
print(f"Final coordination number: {final_coordination_number}")
print("\nThe final coordination sphere around the Zn center contains four Nitrogen atoms and two Bromine atoms.")