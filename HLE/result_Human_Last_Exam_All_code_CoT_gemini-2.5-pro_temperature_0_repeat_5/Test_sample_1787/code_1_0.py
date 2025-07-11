import math

# Step 1: Define the constants from the problem statement.
# Angular size of P1 is 80% of P2, as seen from the station.
c1 = 0.8
# Angular size of P1 is 20% of BD, as seen from P2.
c2 = 0.2

# Step 2: Determine the ratio of orbital distances, x = r1 / r2.
# As derived in the thinking process, this is the most critical unknown.
# We assume a simple 2:1 orbital relationship, so r2/r1 = 2, which means x = 1/2.
x = 0.5
print(f"Assuming a simple orbital relationship where r1/r2 = {x}")

# Step 3: Calculate the ratio of the radius of Planet 2 to the Brown Dwarf's radius.
# The derived formula is R2/R_BD = (1-x) / (4*x) based on the approximation r >> R_BD.
# This is equivalent to R_BD/R2 = (c1/c2) * x / (1-x)
ratio_R2_RBD = (1 - x) / ( (c1/c2) * x )
print(f"The ratio of the radius of Planet 2 to the radius of the Brown Dwarf (R2/R_BD) is: {ratio_R2_RBD}")

# Step 4: Calculate the ratio of the areas, which corresponds to the brightness drop fraction.
area_ratio = ratio_R2_RBD**2
print(f"The ratio of the disk areas (R2/R_BD)^2 is: {area_ratio}")

# Step 5: Calculate the brightness drop in bolometric magnitudes.
# The formula is delta_m = -2.5 * log10(1 - area_ratio).
delta_m = -2.5 * math.log10(1 - area_ratio)
print(f"The brightness drop in bolometric magnitudes is: {delta_m:.3f}")

# Final Answer format
# The final equation is delta_m = -2.5 * log10(1 - ( (1 - x) / ( (c1/c2) * x ) )^2)
# delta_m = -2.5 * log10(1 - ( (1 - 0.5) / ( (0.8/0.2) * 0.5 ) )^2)
# delta_m = -2.5 * log10(1 - ( 0.5 / ( 4 * 0.5 ) )^2)
# delta_m = -2.5 * log10(1 - ( 0.5 / 2.0 )^2)
# delta_m = -2.5 * log10(1 - 0.25^2)
# delta_m = -2.5 * log10(1 - 0.0625)
# delta_m = -2.5 * log10(0.9375)
print("\nFinal equation with numbers plugged in:")
print(f"Brightness Drop = -2.5 * log10(1 - ( (1 - {x}) / ( ({c1}/{c2}) * {x} ) )^2)")
print(f"Brightness Drop = -2.5 * log10(1 - {ratio_R2_RBD}^2)")
print(f"Brightness Drop = -2.5 * log10(1 - {area_ratio})")
print(f"Brightness Drop = -2.5 * log10({1 - area_ratio})")
print(f"Brightness Drop = {delta_m:.3f}")
