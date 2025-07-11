import math

# Define the point of interest P
px, py, pz = 1, -1, 0

# --- Wire 1: Current I along the positive x-axis ---

# The perpendicular distance r1 from P(1, -1, 0) to the x-axis is |-1|
r1 = abs(py)
# The magnitude of the magnetic field B1 is proportional to 1/r1.
# B1_mag = (mu_0 * I) / (2 * pi * r1)
# By the right-hand rule, the direction at P is the -z direction.
# Let's define a base unit B_unit = (mu_0 * I) / (2 * pi).
# B1_vector = - (1/r1) * B_unit in the z-direction.
B1_z_component = -1 / r1

print(f"Analysis for Wire 1 (on x-axis):")
print(f"Point P is at ({px}, {py}, {pz}).")
print(f"Perpendicular distance from P to the x-axis, r1 = {r1}")
print(f"The magnetic field B1 has a component in the z-direction of {B1_z_component:.2f} * (μ₀ * I / (2 * π)).\n")


# --- Wire 2: Current I along the positive y-axis ---

# The perpendicular distance r2 from P(1, -1, 0) to the y-axis is |1|
r2 = abs(px)
# The magnitude of the magnetic field B2 is proportional to 1/r2.
# B2_mag = (mu_0 * I) / (2 * pi * r2)
# By the right-hand rule, the direction at P is the -z direction.
# B2_vector = - (1/r2) * B_unit in the z-direction.
B2_z_component = -1 / r2

print(f"Analysis for Wire 2 (on y-axis):")
print(f"Point P is at ({px}, {py}, {pz}).")
print(f"Perpendicular distance from P to the y-axis, r2 = {r2}")
print(f"The magnetic field B2 has a component in the z-direction of {B2_z_component:.2f} * (μ₀ * I / (2 * π)).\n")


# --- Total Magnetic Field ---

# The total magnetic field is the vector sum. Since both are in the z-direction, we add the components.
B_total_z_component = B1_z_component + B2_z_component
# B_total_z_component is the coefficient for (μ₀ * I / (2 * π))

# The magnitude is the absolute value of this component.
B_total_mag_coeff_2pi = abs(B_total_z_component)
# This means the total magnitude is B_total_mag_coeff_2pi * (μ₀ * I / (2 * π))
# To express it in terms of (μ₀ * I / π), we multiply the coefficient by 2.
B_total_mag_coeff_pi = B_total_mag_coeff_2pi / 2

print("--- Final Result ---")
print(f"The total magnetic field is the sum of the fields from both wires.")
print(f"The total z-component is ({B1_z_component:.2f} + {B2_z_component:.2f}) = {B_total_z_component:.2f} * (μ₀ * I / (2 * π)).")
print(f"This is equivalent to {B_total_z_component/2.0:.2f} * (μ₀ * I / π).")
print("\nThe magnitude is the absolute value of the field.")
print(f"Final Magnitude = {int(abs(B_total_z_component/2.0))} * (μ₀ * I) / ({int(math.pi/math.pi)} * π)")
print("\nThe final equation for the magnitude |B| is:")
print("|B| = (μ₀ * I) / π")
