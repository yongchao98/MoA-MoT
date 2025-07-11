# The genus of the initial surface Sigma.
g_sigma = 10

# We need to find the smallest integer g that serves as the genus for a
# closed surface Sigma' which can contain any surface Sigma with the given properties.
# Sigma' is formed by gluing a "capping surface", S_cap, to the boundary of Sigma.
# The genus of the resulting surface g is the sum of the genus of Sigma and
# the genus of S_cap.
# g = g_sigma + g_cap

# The problem requires this to work for *any* possible embedding of Sigma.
# We must therefore consider the "worst-case" scenario, which maximizes the
# required genus of the capping surface, g_cap.
#
# A surface of genus 10 has 10 handles. The worst-case embedding is one where
# all 10 of these handles are topologically linked with the boundary of Sigma.
# For example, if the boundary is a circle in a plane, each handle can be
# constructed to pass through the disk spanning that circle.
#
# To "cap" such a tangled surface without intersecting it, the capping surface S_cap
# must have a hole for each of Sigma's handles to pass through. A surface with
# N holes (handles) has genus N. Thus, S_cap must have a genus of 10.
g_cap_max = 10

# Calculate the final genus g for this worst-case scenario.
g_final = g_sigma + g_cap_max

print("The genus of the original surface is g_sigma.")
print(f"g_sigma = {g_sigma}")
print("The maximum required genus for the capping surface is g_cap.")
print(f"g_cap = {g_cap_max}")
print("The resulting genus g is the sum of these two values.")
print(f"g = g_sigma + g_cap")
print(f"g = {g_sigma} + {g_cap_max} = {g_final}")
