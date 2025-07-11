import math

# Define the given variables based on the problem description.
# I is the horizontal distance from the gun to the highest point of elevation.
I = 500  # in meters

# According to the principle of the center of mass, the center of mass
# of the fragments will land where the original projectile would have landed.
# Since the trajectory is symmetric, the full range would have been 2 * I.
# This is the landing position of the center of mass (x_cm).
x_cm = 2 * I

# One fragment is stated to have fallen near the gun. We take its landing
# position (x1) to be 0 m.
x1 = 0

# For two fragments of equal mass, the center of mass is the average of their
# landing positions: x_cm = (x1 + x2) / 2.
# We need to find the landing position of the second fragment, x2.
# We can rearrange the formula to solve for x2: x2 = 2 * x_cm - x1.
# Substituting x_cm = 2 * I, we get: x2 = 2 * (2 * I) - x1, which simplifies to x2 = 4 * I.

# Calculate the final maximum distance for the second fragment.
max_distance = 4 * I

# As requested, we print the final equation with the numerical values.
# The formula for the maximum distance (R_max) of the second fragment is R_max = 4 * I.
print(f"The equation for the maximum distance of the second fragment (R_max) is:")
print(f"R_max = 4 * I")
print(f"Substituting I = {I} m:")
print(f"R_max = 4 * {I}")
print(f"R_max = {max_distance} m")
print("\nSo, the maximum distance from the gun where the second fragment can land is 2000 meters.")
