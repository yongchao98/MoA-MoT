import math

# This script calculates the magnitude of the magnetic field at a specific point
# due to two infinite current-carrying wires.

print("### Step-by-step Calculation ###\n")

# 1. Define the point of interest P and constants.
# The problem is symbolic, so we will use mu_0 and I as variables.
x, y, z = 1, -1, 0
print(f"The magnetic field will be calculated at the point P(x,y,z) = ({x}, {y}, {z}).\n")

# 2. Analyze Wire 1 (along x-axis, current in +x direction)
print("--- Wire 1 on x-axis ---")
# The perpendicular distance 'r' from the x-axis to a point (x,y,z) is sqrt(y^2 + z^2).
r1_val = math.sqrt(y**2 + z**2)
print(f"The perpendicular distance from P to Wire 1 is r1 = sqrt(({y})^2 + ({z})^2) = {r1_val}")
# The magnitude of the magnetic field is B = (mu_0 * I) / (2 * pi * r).
# Using the right-hand rule (current in +x), the field at P (y=-1) points into the page (-z direction).
# So, B1_vector = (0, 0, - (mu_0 * I) / (2 * pi * r1))
print(f"The magnetic field B1 has a magnitude of (mu_0 * I) / (2 * pi * {int(r1_val)}).")
print("By the right-hand rule, its direction is along the negative z-axis.")
print("B1_vector = (0, 0, - mu_0*I / (2*pi))\n")


# 3. Analyze Wire 2 (along y-axis, current in +y direction)
print("--- Wire 2 on y-axis ---")
# The perpendicular distance 'r' from the y-axis to a point (x,y,z) is sqrt(x^2 + z^2).
r2_val = math.sqrt(x**2 + z**2)
print(f"The perpendicular distance from P to Wire 2 is r2 = sqrt(({x})^2 + ({z})^2) = {r2_val}")
# The magnitude of the magnetic field is B = (mu_0 * I) / (2 * pi * r).
# Using the right-hand rule (current in +y), the field at P (x=1) also points into the page (-z direction).
# So, B2_vector = (0, 0, - (mu_0 * I) / (2 * pi * r2))
print(f"The magnetic field B2 has a magnitude of (mu_0 * I) / (2 * pi * {int(r2_val)}).")
print("By the right-hand rule, its direction is also along the negative z-axis.")
print("B2_vector = (0, 0, - mu_0*I / (2*pi))\n")


# 4. Calculate the total magnetic field and its magnitude
print("--- Total Magnetic Field ---")
# The total field is the vector sum B_total = B1 + B2.
# Since both vectors point in the -z direction, we add their magnitudes.
# Total z-component = -(mu_0*I / (2*pi)) - (mu_0*I / (2*pi)) = -2*mu_0*I / (2*pi) = -mu_0*I / pi
print("The total field B_total is the vector sum B1 + B2.")
print("B_total = (0, 0, - (mu_0*I / (2*pi)) - (mu_0*I / (2*pi)))")
print("B_total = (0, 0, - (2*mu_0*I) / (2*pi))")
print("B_total = (0, 0, - (mu_0*I) / pi)\n")

# The magnitude is the length of this vector.
print("--- Final Magnitude ---")
print("The magnitude of the total magnetic field is the absolute value of the total vector's length.")
print("Magnitude = | - (mu_0*I) / pi |")

# Final expression with numbers as requested
final_numerator_coefficient = 1
final_denominator_coefficient = 1
print("\nThe final simplified expression for the magnitude is:")
print(f"B = ({final_numerator_coefficient} * mu_0 * I) / ({final_denominator_coefficient} * pi)")
