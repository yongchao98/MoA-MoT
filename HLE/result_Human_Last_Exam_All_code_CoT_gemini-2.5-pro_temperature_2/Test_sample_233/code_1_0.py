# Step 1: Define the genus of the initial surface Σ.
g_sigma = 10

# Step 2: In the worst-case configuration, the number of handles of Σ that
# are topologically linked with its boundary is equal to its genus.
num_linking_handles = g_sigma

# Step 3: The capping surface, C, must have a genus at least equal to the
# number of linking handles to avoid intersecting Σ. This determines the
# required genus of the cap, g(C), in the worst-case scenario.
g_cap_worst_case = num_linking_handles

# Step 4: The genus of the final closed surface Σ' is the sum of the
# genus of the original surface Σ and the genus of the cap C.
# g = g(Σ) + g(C)
g_final = g_sigma + g_cap_worst_case

# Step 5: Print the final calculation, showing each number in the equation.
print(f"The genus of the initial surface Σ is g(Σ) = {g_sigma}.")
print(f"In the worst-case embedding, the capping surface C must have a genus of g(C) = {g_cap_worst_case}.")
print(f"The final genus g is the sum of the genus of Σ and the genus of C.")
print(f"g = g(Σ) + g(C) = {g_sigma} + {g_cap_worst_case} = {g_final}")