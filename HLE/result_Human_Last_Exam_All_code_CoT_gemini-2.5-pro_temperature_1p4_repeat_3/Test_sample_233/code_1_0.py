# The genus of the initial smoothly embedded oriented surface, denoted as Sigma.
g_sigma = 10

# The surface Sigma has a single boundary component which is unknotted.
# In R^3, an unknotted boundary component bounds a smoothly embedded disk.
# A disk is a surface of genus 0. We can use this disk to "cap" Sigma.
g_cap = 0

# To find the genus of the new closed surface (Sigma'), we add the genus of the
# original surface (g_sigma) to the genus of the cap (g_cap).
# Let g be the genus of the smallest such containing surface.
g = g_sigma + g_cap

# The problem asks for the smallest positive integer g such that,
# regardless of our choice of Sigma, there exists a smoothly embedded
# oriented closed surface Sigma' of genus g where Sigma is a subset of Sigma'.
# Since the genus of a surface cannot be less than the genus of its subsurface,
# g must be greater than or equal to g_sigma.
# Our construction shows that g=10 is always possible.
# Thus, the smallest possible value for g is 10.

print("Let g be the genus of the final closed surface Sigma'.")
print("Let g_sigma be the genus of the initial surface Sigma.")
print("Let g_cap be the genus of the surface used to cap the boundary.")
print("\nThe final genus g is the sum of the initial genus and the cap's genus:")
print(f"g = g_sigma + g_cap")
print(f"g = {g_sigma} + {g_cap}")
print(f"g = {g}")
print("\nThe smallest positive integer g is therefore 10.")