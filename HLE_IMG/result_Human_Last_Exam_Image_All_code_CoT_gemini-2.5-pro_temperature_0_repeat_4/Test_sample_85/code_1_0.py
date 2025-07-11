import math

# This script outlines the calculation for the furthest distance on the cone's surface.
# The problem is solved by unrolling the cone into a semicircle of radius 'd'.
# The starting point P is on the arc of the semicircle, and the furthest point Q
# is at the end of the semicircle's diameter.
# We use the Pythagorean theorem to find the distance in the unrolled 2D plane.

print("The equation for the square of the distance is:")
print("Distance^2 = (d - 0)^2 + (0 - d)^2")
print("Distance^2 = d^2 + d^2")

# Define the coefficient for d^2
coefficient = 2
print(f"Distance^2 = {coefficient} * d^2")

print("\nTaking the square root gives the final distance:")
# The number in the final equation is 2
final_number = 2
print(f"Distance = d * sqrt({final_number})")