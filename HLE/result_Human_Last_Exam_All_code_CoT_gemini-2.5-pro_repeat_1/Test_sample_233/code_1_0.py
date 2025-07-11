# The genus of the initial surface Sigma.
genus_sigma = 10

# In the worst-case scenario, each of the handles of Sigma
# can be tangled in a way that requires a corresponding handle
# on the capping surface. The number of handles is equal to the genus.
# Therefore, the genus of the capping surface must be 10.
genus_cap_worst_case = 10

# The genus of the final closed surface Sigma' is the sum of the
# genus of the original surface and the genus of the capping surface.
final_genus = genus_sigma + genus_cap_worst_case

# Print the equation and the result.
print(f"The genus of the original surface Σ is {genus_sigma}.")
print(f"The required genus of the capping surface in the worst case is {genus_cap_worst_case}.")
print(f"The smallest positive integer g for the final surface Σ' is the sum of these two genera.")
print(f"g = {genus_sigma} + {genus_cap_worst_case} = {final_genus}")