import math

# This script calculates the vertical velocity of the drawbridge edge.

# Step 1: Define the constants and values needed for the calculation.
# The derived formula for the vertical velocity (dh/dt) is:
# dh/dt = (3 * pi * sqrt(15)) / cos(pi/12)

val_3 = 3
val_pi = math.pi
val_15 = 15
val_12 = 12

# Step 2: Calculate the value of each term in the formula.
sqrt_15_val = math.sqrt(val_15)
cos_pi_div_12_val = math.cos(val_pi / val_12)

# Step 3: Print the final equation showing each numerical value.
# The prompt asks to show each number in the final equation.
print("The final equation for the vertical velocity (dh/dt) is:")
print("dh/dt = (3 * π * √15) / cos(π/12)")
print("\nPlugging in the computed values for the terms:")
print(f"dh/dt = ({val_3} * {val_pi:.5f} * {sqrt_15_val:.5f}) / {cos_pi_div_12_val:.5f}")

# Step 4: Calculate and print the final result.
# Assuming units are m/min based on context in the problem description.
vertical_velocity = (val_3 * val_pi * sqrt_15_val) / cos_pi_div_12_val
print(f"\nThe calculated vertical velocity is {vertical_velocity:.4f} m/min.")

# Also show the velocity in m/s for reference.
velocity_mps = vertical_velocity / 60
print(f"This is equivalent to {velocity_mps:.4f} m/s.")