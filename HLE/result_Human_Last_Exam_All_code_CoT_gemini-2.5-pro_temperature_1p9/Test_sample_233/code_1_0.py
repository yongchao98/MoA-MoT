# The genus of the given surface Sigma.
genus_sigma = 10

# In the worst-case scenario, the surface Sigma can be embedded in R^3
# such that all 10 of its handles are topologically linked with its boundary.
# To create a closed surface Sigma' that contains Sigma, we must cap the boundary
# with another surface, F.
# For F to not intersect the 10 linking handles of Sigma, it must have
# a topological complexity (genus) of its own. It requires one handle for each
# handle of Sigma it needs to avoid.
genus_f_worst_case = 10

# The genus of the final surface, g, is the sum of the genus of the original
# surface and the genus of the capping surface.
g = genus_sigma + genus_f_worst_case

print(f"The genus of the initial surface Σ is {genus_sigma}.")
print(f"The necessary genus of the capping surface F in the worst-case scenario is {genus_f_worst_case}.")
print(f"The smallest positive integer g such that Σ' can contain any such Σ is:")
print(f"g = {genus_sigma} + {genus_f_worst_case} = {g}")