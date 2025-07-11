import math

# This script calculates the furthest distance from a point A to the surface of a
# 1 cubic meter object shaped to maximize the gravitational field at A.

# The final equation for this distance (r_max) is derived from physics principles:
# r_max = (numerator / (denominator_const * pi)) ^ exponent

# --- Numbers in the final equation ---
# This is derived from the volume calculation of the optimal shape.
numerator = 15
# This also comes from the volume calculation.
denominator_const = 4
# The mathematical constant pi.
pi = math.pi
# The exponent is 1/3, resulting from solving for the maximum radius.
exponent = 1/3

# --- Calculation ---
# We compute the value based on the equation r_max = (15 / (4 * pi))^(1/3).
r_max = (numerator / (denominator_const * pi)) ** exponent

# --- Output ---
# The following lines print the equation with its components and the final result.
print("To find the furthest point, we first derive its formula.")
print(f"The furthest distance, r_max, is calculated using the equation:")
print(f"r_max = (numerator / (denominator_const * pi)) ^ exponent")
print(f"where:")
print(f"  numerator = {numerator}")
print(f"  denominator_const = {denominator_const}")
print(f"  pi = {pi}")
print(f"  exponent = {exponent}")
print("\nCalculating the result:")
print(f"The furthest point on the surface of the playdough is {r_max:.6f} meters from point A.")
