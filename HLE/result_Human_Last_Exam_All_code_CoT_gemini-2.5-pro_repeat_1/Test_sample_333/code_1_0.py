import math

# Horizontal distance from the gun to the highest point of elevation in meters.
I = 500

# The first fragment fell near the gun, so its landing distance is 0.
x1 = 0

# The center of mass (CM) of the fragments continues on the original trajectory.
# The total range of this trajectory (the landing spot of the CM) is twice
# the distance to the highest point.
R_cm = 2 * I

# For two fragments of equal mass, the landing position of the CM is the
# average of their landing positions: R_cm = (x1 + x2) / 2.
# We solve for x2, the landing position of the second fragment: x2 = 2 * R_cm - x1.
x2 = 2 * R_cm - x1

# The problem asks to output the final equation with each number.
# The equation is: 2 * R_cm - x1 = x2
print(f"The maximum safe distance can be calculated from the equation: 2 * (Range of Center of Mass) - (Landing position of first fragment) = (Landing position of second fragment)")
print(f"The final calculation is:")
print(f"{2} * {R_cm} - {x1} = {int(x2)}")
print(f"\nTherefore, the maximum safe distance from the gun is {int(x2)} meters.")
