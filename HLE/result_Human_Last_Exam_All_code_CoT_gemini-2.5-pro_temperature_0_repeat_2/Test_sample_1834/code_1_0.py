import math

# This script calculates the magnitude of the magnetic field from two perpendicular wires.

# Plan:
# 1. Define the coordinates of the point P.
# 2. Calculate the magnetic field contribution from Wire 1 (on the x-axis).
# 3. Calculate the magnetic field contribution from Wire 2 (on the y-axis).
# 4. Sum the contributions and find the magnitude of the total magnetic field.
# 5. Print the final equation, showing the numbers used in the calculation.

# The magnitude of the magnetic field from an infinite wire is B = (mu_0 * I) / (2 * pi * r)
# where r is the perpendicular distance from the wire.

# 1. Define the point P
x, y, z = 1, -1, 0

# 2. Contribution from Wire 1 (on x-axis, current in +x)
# The perpendicular distance from the x-axis to P(x,y,z) is r1 = sqrt(y^2 + z^2)
r1 = math.sqrt(y**2 + z**2)
# By the right-hand rule, the field is in the -z direction.
# The magnitude is B1 = (mu_0 * I) / (2 * pi * r1)

# 3. Contribution from Wire 2 (on y-axis, current in +y)
# The perpendicular distance from the y-axis to P(x,y,z) is r2 = sqrt(x^2 + z^2)
r2 = math.sqrt(x**2 + z**2)
# By the right-hand rule, the field is also in the -z direction.
# The magnitude is B2 = (mu_0 * I) / (2 * pi * r2)

# 4. Total magnetic field and its magnitude
# Since both B1 and B2 point in the same direction (-z), the total magnitude is the sum of their individual magnitudes.
# B_total = B1 + B2

# 5. Print the final equation with numbers
print("The magnitude of the total magnetic field B_total is the sum of the magnitudes from each wire.")

print("\nContribution from Wire 1:")
print(f"The distance r1 from the x-axis to ({x}, {y}, {z}) is {r1}.")
print(f"The magnitude B1 = (mu_0 * I) / (2 * pi * r1) = (mu_0 * I) / (2 * pi * {r1})")

print("\nContribution from Wire 2:")
print(f"The distance r2 from the y-axis to ({x}, {y}, {z}) is {r2}.")
print(f"The magnitude B2 = (mu_0 * I) / (2 * pi * r2) = (mu_0 * I) / (2 * pi * {r2})")

print("\nTotal Magnitude:")
print("B_total = B1 + B2")
# The numbers in the equation are 2, pi, r1, and r2.
print(f"B_total = (mu_0 * I) / (2 * pi * {r1}) + (mu_0 * I) / (2 * pi * {r2})")
print(f"B_total = 2 * (mu_0 * I) / (2 * pi * {r1})") # Since r1=r2=1
print(f"B_total = (mu_0 * I) / (pi * {r1})")

print("\nSince r1 = 1, the final equation for the magnitude is:")
print("B_total = (mu_0 * I) / pi")