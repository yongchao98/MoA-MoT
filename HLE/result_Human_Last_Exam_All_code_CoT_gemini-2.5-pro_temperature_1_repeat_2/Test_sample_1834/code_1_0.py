import math

# The problem is to find the magnitude of the magnetic field at point P(1, -1, 0)
# due to two infinite wires.
# Wire 1: along x-axis, current I in +x direction.
# Wire 2: along y-axis, current I in +y direction.

# The magnetic field B from a long straight wire is given by B = (mu_0 * I) / (2 * pi * r),
# where r is the perpendicular distance from the wire.
# The total field is the vector sum of the fields from each wire: B_total = B1 + B2.

# Define the point of interest
x, y, z = 1, -1, 0

# We will express the magnetic field in terms of the symbolic unit B_unit = (μ₀ * I) / (2π).
# The coefficients of this unit will be calculated.

# --- Field from Wire 1 (on x-axis) ---

# 1. Calculate the perpendicular distance r1 from P(x, y, z) to the x-axis.
r1 = math.sqrt(y**2 + z**2)

# 2. The magnitude B1 is proportional to 1/r1.
B1_mag_coeff = 1.0 / r1

# 3. Determine the direction of B1 using the right-hand rule.
# Current is in +x. Point P is at y=-1, z=0.
# The field lines are circles around the x-axis. At P, the field points in the -z direction.
# So, B1_vector = (0, 0, -B1_mag_coeff) in units of B_unit.
B1_vec = (0, 0, -B1_mag_coeff)

# --- Field from Wire 2 (on y-axis) ---

# 1. Calculate the perpendicular distance r2 from P(x, y, z) to the y-axis.
r2 = math.sqrt(x**2 + z**2)

# 2. The magnitude B2 is proportional to 1/r2.
B2_mag_coeff = 1.0 / r2

# 3. Determine the direction of B2 using the right-hand rule.
# Current is in +y. Point P is at x=1, z=0.
# The field lines are circles around the y-axis. At P, the field points in the +z direction.
# So, B2_vector = (0, 0, +B2_mag_coeff) in units of B_unit.
B2_vec = (0, 0, B2_mag_coeff)

# --- Total Magnetic Field ---

# Calculate the total magnetic field vector by summing the individual vectors.
B_total_vec_x = B1_vec[0] + B2_vec[0]
B_total_vec_y = B1_vec[1] + B2_vec[1]
B_total_vec_z = B1_vec[2] + B2_vec[2]

# Calculate the magnitude of the total magnetic field.
# The final magnitude is sqrt(Bx^2 + By^2 + Bz^2) * B_unit
total_mag_coeff = math.sqrt(B_total_vec_x**2 + B_total_vec_y**2 + B_total_vec_z**2)

# --- Print the Calculation and Result ---
print("This script calculates the magnitude of the magnetic field at P(1, -1, 0).")
print("The result is expressed in terms of the constant B_unit = (μ₀ * I) / (2π).\n")

print(f"Field from Wire 1 (x-axis):")
print(f"  - Distance r1 = sqrt({y}^2 + {z}^2) = {r1}")
print(f"  - B1 Vector = (0, 0, -1/{r1}) * B_unit = (0, 0, {B1_vec[2]}) * B_unit")

print(f"\nField from Wire 2 (y-axis):")
print(f"  - Distance r2 = sqrt({x}^2 + {z}^2) = {r2}")
print(f"  - B2 Vector = (0, 0, +1/{r2}) * B_unit = (0, 0, {B2_vec[2]}) * B_unit")

print("\nTotal Magnetic Field Vector B_total = B1 + B2:")
print(f"  - B_total = (0, 0, {B1_vec[2]} + {B2_vec[2]}) * B_unit = ({B_total_vec_x}, {B_total_vec_y}, {B_total_vec_z}) * B_unit")

print("\nMagnitude of the Total Magnetic Field |B_total|:")
print(f"  - |B_total| = sqrt({B_total_vec_x}^2 + {B_total_vec_y}^2 + {B_total_vec_z}^2) * B_unit")
print(f"  - |B_total| = {total_mag_coeff} * B_unit")

print(f"\nThe final magnitude of the magnetic field is: {total_mag_coeff}")