# The genus of the given surface Sigma.
genus_sigma = 10

# The problem asks for the smallest genus g of a closed surface Sigma'
# such that Sigma is a subsurface of Sigma'.
# The relationship between the genera is g(Sigma') = g(Sigma) + g(F),
# where F is the "capping" surface.
# We need to find the maximum possible minimal genus for F, over all
# possible embeddings of Sigma.
# Let g_cap_max be this value.
# Based on results in 3D topology, the maximum required genus for the
# capping surface is equal to the genus of the original surface.
g_cap_max = genus_sigma

# The final genus g is the sum of the genus of Sigma and this maximum
# capping genus.
g = genus_sigma + g_cap_max

# Output the final equation
print(f"The genus of the initial surface Σ is g(Σ) = {genus_sigma}.")
print(f"The maximum required genus for the capping surface F is g_cap_max = {g_cap_max}.")
print(f"The resulting genus g is calculated by the formula g = g(Σ) + g_cap_max.")
print(f"So, g = {genus_sigma} + {g_cap_max} = {g}")