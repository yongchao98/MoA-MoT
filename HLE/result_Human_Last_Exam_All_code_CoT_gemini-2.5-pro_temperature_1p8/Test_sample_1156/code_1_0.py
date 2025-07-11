import math

# Calculate the Golden Ratio G
G = (1 + math.sqrt(5)) / 2

# Calculate the normalization constant C
# I = 2 * (1 - G * ln(G))
# C = 1 / I
C = 1 / (2 * (1 - G * math.log(G)))

# Print the final equation for the normalized density rho(x)
print("The normalised density of the invariant measure is rho(x) = C / (G + sqrt(x))")
print("where:")
print(f"G (Golden Ratio) = {G}")
print(f"C (Normalization Constant) = {C}")
print("\nSo the final equation is:")
print(f"rho(x) = {C} / ({G} + sqrt(x))")
