import math

# --- Problem Description ---
# We are calculating the furthest point on the surface of a 1 cubic meter volume
# that has been shaped to create the strongest possible gravitational field at a point A.
# Let this maximum distance be r_max.

# --- Formula Derivation ---
# The optimal shape has a volume V related to its maximum dimension r_max by the formula:
# V = (4 * pi * r_max^3) / 15
# We can rearrange this to solve for r_max:
# r_max = (15 * V / (4 * pi))^(1/3)

# --- Calculation ---

# Define the given and constant values
volume_V = 1.0
constant_15 = 15.0
constant_4 = 4.0
pi = math.pi
exponent = 1/3

# Perform the calculation
numerator = constant_15 * volume_V
denominator = constant_4 * pi
r_max_cubed = numerator / denominator
r_max = r_max_cubed**exponent

# --- Output ---
# As requested, we will output each number in the final equation.
print("The equation to find the furthest distance (r_max) is:")
print("r_max = (15 * V / (4 * pi))^(1/3)")
print("\nWhere the values are:")
print(f"V = {volume_V}")
print(f"15 = {constant_15}")
print(f"4 = {constant_4}")
print(f"pi = {pi}")
print(f"The exponent (1/3) is {exponent}")
print("\nSubstituting these values into the equation:")
print(f"r_max = ({constant_15} * {volume_V} / ({constant_4} * {pi})) ^ ({exponent})")

# Print the final answer
print(f"\nThe furthest point on the surface of the playdough is {r_max} meters from point A.")
