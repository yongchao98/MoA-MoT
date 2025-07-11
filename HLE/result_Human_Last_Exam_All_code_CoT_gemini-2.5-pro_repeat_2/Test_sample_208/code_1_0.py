import math

# This script calculates the radius of a circle that tightly packs 14 unit circles.
# The formula for this specific packing problem (n=14) is R = 3 + sec(pi/9).

print("To find the radius of a circle that tightly packs fourteen circles of radius one, we use a known mathematical formula.")
print("For 14 circles, the radius 'R' of the enclosing circle is given by R = 3 + sec(pi/9).")
print("This is equivalent to R = 3 + (1 / cos(pi/9)).")

# Define the numbers in the final equation
val_3 = 3
val_1 = 1
angle_rad = math.pi / 9
cos_val = math.cos(angle_rad)

# The final radius R
R = val_3 + (val_1 / cos_val)

print("\nHere are the individual numbers in the final equation:")
print(f"The number three: {val_3}")
print(f"The number one: {val_1}")
# We use the calculated value for cos(pi/9) as a number in the equation.
print(f"The value of cos(pi/9): {cos_val:.8f}")

print("\nPerforming the calculation step-by-step:")
print(f"R = {val_3} + (1 / {cos_val:.8f})")
print(f"R = {val_3} + {1/cos_val:.8f}")
print(f"R = {R:.8f}")

# The '.4g' format specifier rounds the number to 4 significant digits.
print(f"\nThe radius of this circle up to 4 significant digits is {R:.4g}")