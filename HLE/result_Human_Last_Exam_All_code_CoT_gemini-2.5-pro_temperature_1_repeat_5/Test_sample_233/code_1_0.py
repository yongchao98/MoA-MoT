# The genus of the given surface Sigma.
g_sigma = 10

# The problem states Sigma has a single boundary component.
b_sigma = 1

# The problem asks for the smallest genus g for a closed surface Sigma'
# that contains Sigma, regardless of how Sigma is embedded.
# This surface Sigma' is formed by attaching a capping surface F to the boundary of Sigma.
# The genus of the resulting surface is g(Sigma') = g(Sigma) + g(F).

# The challenge is that Sigma can be embedded in a "knotted" way.
# This knottedness can force the capping surface F to have a non-zero genus
# to avoid intersecting Sigma.
# The maximum necessary genus for the capping surface F is bounded by the genus of Sigma itself.
# In the worst-case scenario, we need one handle on F for each handle on Sigma.
g_F_max = g_sigma

# Therefore, the smallest g that works for all possible embeddings is the one
# calculated from this worst-case scenario.
g_final = g_sigma + g_F_max

print("The genus of the final surface, g, is the sum of the genus of the initial surface, g(Sigma), and the maximum necessary genus of the capping surface, g_max(F).")
print(f"g(Sigma) = {g_sigma}")
print(f"g_max(F) = {g_F_max}")
print("The final equation is:")
print(f"{g_sigma} + {g_F_max} = {g_final}")
