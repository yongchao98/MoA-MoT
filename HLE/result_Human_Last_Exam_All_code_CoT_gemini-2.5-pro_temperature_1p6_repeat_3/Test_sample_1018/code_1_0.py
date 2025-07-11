import math

# Define the genera of the surfaces
g1 = 31
g2 = 17

# Dimensions of the surfaces (both are 2-dimensional)
m = 2
n = 2

# Formula for the simplicial volume of a surface Sigma_g for g >= 2
def simplicial_volume_surface(g):
  if g < 2:
    return 0
  return 4 * g - 4

# Calculate the simplicial volume for each surface
vol1 = simplicial_volume_surface(g1)
vol2 = simplicial_volume_surface(g2)

# Calculate the binomial coefficient C(m+n, m)
total_dim = m + n
binom_coeff = math.comb(total_dim, m)

# Apply Gromov's product formula
total_volume = binom_coeff * vol1 * vol2

# Print the final equation with all the computed values
print("The simplicial volume of Σ_{31} x Σ_{17} is calculated using Gromov's product formula:")
print("||Σ_{31} x Σ_{17}|| = C(dim(Σ_{31}) + dim(Σ_{17}), dim(Σ_{31})) * ||Σ_{31}|| * ||Σ_{17}||")
print(f"\nWhere:")
print(f"||Σ_{31}|| = 4*{g1} - 4 = {vol1}")
print(f"||Σ_{17}|| = 4*{g2} - 4 = {vol2}")
print(f"C({total_dim}, {m}) = {binom_coeff}")
print("\nFinal calculation:")
print(f"{binom_coeff} * {vol1} * {vol2} = {total_volume}")