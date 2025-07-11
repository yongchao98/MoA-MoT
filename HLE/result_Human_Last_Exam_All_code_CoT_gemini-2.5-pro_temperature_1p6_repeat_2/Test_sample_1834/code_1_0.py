import math

# The problem asks for a symbolic answer, so we will represent the constants
# mu_0 (permeability of free space) and I (current) as strings for the final output.
mu_0 = "μ₀"
I = "I"
pi = "π"

# Define the point P where the magnetic field is to be calculated.
P = (1, -1, 0)
x, y, z = P

print("Calculating the magnetic field at point P(1, -1, 0).\n")

# --- Contribution from Wire 1 (along x-axis, current in +x) ---

# The distance r1 from the x-axis to point P(x,y,z) is sqrt(y^2 + z^2).
r1 = math.sqrt(y**2 + z**2)
print(f"Wire 1 is on the x-axis.")
print(f"The perpendicular distance from Wire 1 to point P is r1 = sqrt(({y})^2 + ({z})^2) = {r1}.")

# The magnitude of the magnetic field from Wire 1 is B1 = (mu_0 * I) / (2 * pi * r1).
# By the right-hand rule, with current in +x, the field at P (with y=-1) points in the -z direction.
B1_direction = "-z"
print(f"The magnetic field B1 from Wire 1 points in the {B1_direction} direction.")
print(f"The magnitude of B1 is |B1| = ({mu_0} * {I}) / (2 * {pi} * {r1}) = ({mu_0}*{I})/(2{pi}).\n")

# --- Contribution from Wire 2 (along y-axis, current in +y) ---

# The distance r2 from the y-axis to point P(x,y,z) is sqrt(x^2 + z^2).
r2 = math.sqrt(x**2 + z**2)
print(f"Wire 2 is on the y-axis.")
print(f"The perpendicular distance from Wire 2 to point P is r2 = sqrt(({x})^2 + ({z})^2) = {r2}.")

# The magnitude of the magnetic field from Wire 2 is B2 = (mu_0 * I) / (2 * pi * r2).
# By the right-hand rule, with current in +y, the field at P (with x=1) points in the -z direction.
B2_direction = "-z"
print(f"The magnetic field B2 from Wire 2 points in the {B2_direction} direction.")
print(f"The magnitude of B2 is |B2| = ({mu_0} * {I}) / (2 * {pi} * {r2}) = ({mu_0}*{I})/(2{pi}).\n")

# --- Total Magnetic Field ---

print("Both B1 and B2 point in the same direction (-z).")
print("Therefore, the total magnetic field magnitude is the sum of their individual magnitudes.\n")

# B_total = B1 + B2
print(f"Final Calculation:")
print(f"|B_total| = |B1| + |B2|")
print(f"|B_total| = ({mu_0} * {I}) / (2 * {pi} * {r1}) + ({mu_0} * {I}) / (2 * {pi} * {r2})")
print(f"|B_total| = 2 * (({mu_0} * {I}) / (2 * {pi}))")
print(f"|B_total| = ({mu_0} * {I}) / {pi}")
