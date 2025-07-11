import math

# Define the point P(x, y, z)
x, y, z = 1, -1, 0

print("This script calculates the magnetic field at a point P(x, y, z) due to two infinite wires.")
print("The magnitude of the magnetic field from a single infinite wire is given by B = (mu_0 * I) / (2 * pi * r).")
print("For simplicity, let's define a constant C = mu_0 * I / (2 * pi). Then, the magnitude is B = C / r, where r is the perpendicular distance to the wire.\n")

print(f"The point of interest is P(x,y,z) = ({x}, {y}, {z}).\n")

# --- Analysis for Wire 1 ---
print("--- Wire 1: Placed along the x-axis, with current in the +x direction ---")
# The perpendicular distance from P to the x-axis is r1.
r1 = math.sqrt(y**2 + z**2)
print(f"The perpendicular distance r1 from point P to Wire 1 is given by sqrt(y^2 + z^2).")
print(f"r1 = sqrt(({y})^2 + ({z})^2) = {r1:.2f}")

print("The magnitude of the magnetic field from this wire, B1, is C / r1.")
print("Using the right-hand rule (with thumb in +x direction), the field at P(1, -1, 0) points in the +z direction.")
# The z-component of B1 will be C / r1. We store the coefficient of C.
b1_z_coeff = 1 / r1
print(f"So, the magnetic field vector from Wire 1 is B1 = (0, 0, C/r1) = (0, 0, C/{r1:.2f}) = (0, 0, {b1_z_coeff:.2f}*C).\n")

# --- Analysis for Wire 2 ---
print("--- Wire 2: Placed along the y-axis, with current in the +y direction ---")
# The perpendicular distance from P to the y-axis is r2.
r2 = math.sqrt(x**2 + z**2)
print(f"The perpendicular distance r2 from point P to Wire 2 is given by sqrt(x^2 + z^2).")
print(f"r2 = sqrt(({x})^2 + ({z})^2) = {r2:.2f}")

print("The magnitude of the magnetic field from this wire, B2, is C / r2.")
print("Using the right-hand rule (with thumb in +y direction), the field at P(1, -1, 0) points in the -z direction.")
# The z-component of B2 will be -C / r2. We store the coefficient of C.
b2_z_coeff = -1 / r2
print(f"So, the magnetic field vector from Wire 2 is B2 = (0, 0, -C/r2) = (0, 0, -C/{r2:.2f}) = (0, 0, {b2_z_coeff:.2f}*C).\n")

# --- Superposition ---
print("--- Total Magnetic Field ---")
print("By the principle of superposition, the total magnetic field, B_total, is the vector sum of B1 and B2.")
print("B_total = B1 + B2")
# The x and y components of both B1 and B2 are zero. We only need to sum the z-components.
total_z_coeff = b1_z_coeff + b2_z_coeff
print(f"The equation for the z-component is:")
print(f"B_total_z = B1_z + B2_z = ({b1_z_coeff:.2f} * C) + ({b2_z_coeff:.2f} * C)")
print(f"B_total_z = ({b1_z_coeff:.2f} + {b2_z_coeff:.2f}) * C = {total_z_coeff:.2f} * C")
print(f"Thus, the total magnetic field vector B_total is (0, 0, 0).\n")

# --- Final Magnitude ---
print("--- Final Magnitude Calculation ---")
# The magnitude of the vector (0, 0, 0) is 0.
final_magnitude = math.sqrt(0**2 + 0**2 + 0**2)
print("The magnitude of the total magnetic field is the magnitude of the B_total vector.")
print(f"|B_total| = sqrt((0*C)^2 + (0*C)^2 + ({total_z_coeff:.2f}*C)^2) = {final_magnitude:.2f}")
print("\nTherefore, the final magnitude of the magnetic field is 0.")