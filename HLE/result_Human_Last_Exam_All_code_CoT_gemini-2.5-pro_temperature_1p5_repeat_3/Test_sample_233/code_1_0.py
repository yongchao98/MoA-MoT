# The problem asks for the smallest integer g that is the genus of a closed surface Sigma'
# which contains an embedded surface Sigma of genus 10 with one unknotted boundary component.
# This value g must work regardless of how Sigma is embedded.

# Genus of the initial surface Sigma
genus_sigma = 10

# The final surface Sigma' is formed by capping the boundary of Sigma with a surface D.
# The genus of the final surface g is the sum of the genus of Sigma and the genus of the cap D.
# g = genus_sigma + genus_d

# We need to find the g that works for the "worst-case" embedding of Sigma.
# A surface of genus 10 has 10 "handles". The worst-case embedding is when these 10 handles
# obstruct the capping surface. Each obstruction requires adding a handle to the cap to bypass it.
# Therefore, the maximum required genus for the capping surface D is equal to the genus of Sigma.
max_genus_d = 10

# The smallest possible value for g that works for all cases is determined by this worst-case scenario.
g = genus_sigma + max_genus_d

# We print the final equation as requested.
print("The genus 'g' of the final surface Sigma' is determined by the equation:")
print(f"g = genus(Sigma) + genus(Capping Surface)")
print(f"{g} = {genus_sigma} + {max_genus_d}")