import math

# Step 1: Derivations from the problem description.
#
# From the first observation (from the space station):
# Angular size of P1 = 0.8 * Angular size of P2
# R_1 / d_1 ≈ 0.8 * R_2 / d_2  (using approximations d >> R_BD)
#
# From the second observation (from Planet 2):
# Angular size of P1 = 0.2 * Angular size of the brown dwarf
# R_1 / (d_2 - d_1) = 0.2 * R_BD / d_2
#
# Let x = d_1 / d_2 and k = R_2 / R_BD. Combining the two observations gives:
# k = (1 - x) / (4 * x)

# Step 2: Use orbital mechanics to find x.
# Velocity of P1 (parabolic orbit, at pericenter d_1): v_1 = sqrt(2GM/d_1)
# Velocity of P2 (circular orbit, at radius d_2): v_2 = sqrt(GM/d_2)
# Assume the simple relation v_1 = 2 * v_2.
# sqrt(2GM/d_1) = 2 * sqrt(GM/d_2)
# Squaring both sides: 2GM/d_1 = 4GM/d_2
# This simplifies to 2/d_1 = 4/d_2, which means d_1/d_2 = 2/4 = 0.5.
# So, x = 0.5.

# Step 3: Calculate k using the value of x.
# k = (1 - 0.5) / (4 * 0.5) = 0.5 / 2.0 = 0.25
# This means the radius of Planet 2 is 1/4 the radius of the brown dwarf.

# Step 4: Calculate the magnitude drop.
# The formula is delta_m = -2.5 * log10(1 - k^2).
# k^2 = (1/4)^2 = 1/16.
# So, delta_m = -2.5 * log10(1 - 1/16) = -2.5 * log10(15/16).

# Final calculation, showing each number in the equation.
constant_factor = 2.5
numerator = 15
denominator = 16

fraction = numerator / denominator
magnitude_drop = -constant_factor * math.log10(fraction)

print(f"The calculation for the magnitude drop is: delta_m = -{constant_factor} * log10({numerator}/{denominator})")
print(f"The brightness drop of the brown dwarf is {magnitude_drop:.3f} magnitudes.")

<<<0.070>>>