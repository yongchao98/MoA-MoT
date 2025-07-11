import math

# The problem is to find the maximum distance from a point A to the surface
# of a 1 cubic meter body that creates the strongest possible gravitational field at A.

# The optimal shape is bounded by the surface r^2 = R^2 * cos(θ).
# The furthest point from A on this surface is at a distance R.
# The volume (V) of this shape is given by the formula: V = (4 * π * R^3) / 15.

# We are given V = 1 cubic meter. We need to solve for R.
# 1 = (4 * π * R^3) / 15
# R^3 = (15 * 1) / (4 * π)
# R = ((15 * 1) / (4 * π))^(1/3)

# Define the constants from the equation.
volume = 1.0
numerator_val = 15.0
denominator_val_1 = 4.0
denominator_val_2 = math.pi

# Calculate R, the maximum distance.
R = (numerator_val * volume / (denominator_val_1 * denominator_val_2))**(1.0/3.0)

# Print the final equation with the numbers plugged in, as requested.
print("The furthest distance (R) is calculated using the volume formula for the optimal shape.")
print("The final equation is: R = (15 * V / (4 * π))^(1/3)")
print("\nPlugging in the given volume V = 1.0 m^3:")
print(f"R = ({numerator_val} * {volume} / ({denominator_val_1} * {denominator_val_2}))^(1/3)")
print("\nResult:")
print(f"The furthest point on the surface of the playdough is {R:.6f} meters from point A.")
