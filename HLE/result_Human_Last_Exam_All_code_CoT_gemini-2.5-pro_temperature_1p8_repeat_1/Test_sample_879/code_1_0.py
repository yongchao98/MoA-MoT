# Define the known masses from the experimental data.
kag1_monomer_mass = 32350
complex_mass_in_og = 101553
protein_subunits_in_complex = 3
lipid_subunits_in_complex = 3

# Calculate the mass of the protein trimer.
protein_trimer_mass = protein_subunits_in_complex * kag1_monomer_mass

# Calculate the total mass of the bound non-protein molecules.
total_lipid_mass = complex_mass_in_og - protein_trimer_mass

# Calculate the mass of a single lipid molecule.
single_lipid_mass = total_lipid_mass / lipid_subunits_in_complex

print("Analysis of the Kag1 complex composition:")
print(f"The calculation demonstrates that the complex is a trimer of Kag1 stabilized by three lipid molecules.")
print("The final equation is:")
print(f"{protein_subunits_in_complex} * {kag1_monomer_mass} Da (Kag1) + {lipid_subunits_in_complex} * {int(round(single_lipid_mass))} Da (Lipid) = {complex_mass_in_og} Da")