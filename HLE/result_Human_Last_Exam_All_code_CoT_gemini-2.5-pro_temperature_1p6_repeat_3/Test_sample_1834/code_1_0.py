import numpy as np

# This script calculates the magnitude of the magnetic field at a specific point
# due to two infinite current-carrying wires.

# --- Problem Setup ---
# Wire 1: along x-axis, current I in +x direction.
# Wire 2: along y-axis, current I in +y direction.
# Point P: (x, y, z) = (1, -1, 0).
# The final answer is expressed in terms of mu_0 and I.

# Define the point of interest
P = np.array([1, -1, 0])
x, y, z = P

print(f"Calculating the magnetic field at point P(x, y, z) = ({x}, {y}, {z}).")
print("-" * 50)

# --- Contribution from Wire 1 (along x-axis) ---
print("Step 1: Calculate the magnetic field from Wire 1 (on x-axis).")

# The formula for the magnetic field magnitude from an infinite wire is B = (mu_0 * I) / (2 * pi * r).
# For Wire 1 on the x-axis, the perpendicular distance to P(x, y, z) is r1 = sqrt(y^2 + z^2).
r1 = np.sqrt(y**2 + z**2)
print(f"The perpendicular distance from P to Wire 1 is r1 = sqrt({y}^2 + {z}^2) = {r1:.2f}")

# The direction is given by the right-hand rule. For a current in +x, at a point with y=-1,
# the magnetic field points in the -z direction.
# Vector B1 will be (0, 0, -B1_magnitude).
print(f"By the right-hand rule, the B-field from Wire 1 at P is in the -z direction.")
print(f"The magnitude of the field is |B1| = mu_0 * I / (2 * pi * {r1:.2f})")
print("-" * 50)


# --- Contribution from Wire 2 (along y-axis) ---
print("Step 2: Calculate the magnetic field from Wire 2 (on y-axis).")

# For Wire 2 on the y-axis, the perpendicular distance to P(x, y, z) is r2 = sqrt(x^2 + z^2).
r2 = np.sqrt(x**2 + z**2)
print(f"The perpendicular distance from P to Wire 2 is r2 = sqrt({x}^2 + {z}^2) = {r2:.2f}")

# The direction is given by the right-hand rule. For a current in +y, at a point with x=1,
# the magnetic field points in the -z direction.
# Vector B2 will be (0, 0, -B2_magnitude).
print(f"By the right-hand rule, the B-field from Wire 2 at P is also in the -z direction.")
print(f"The magnitude of the field is |B2| = mu_0 * I / (2 * pi * {r2:.2f})")
print("-" * 50)

# --- Total Magnetic Field ---
print("Step 3: Calculate the total magnetic field by vector addition.")

# Since both B1 and B2 point in the same direction (-z), we can add their magnitudes.
print("B_total = B1 + B2")
print("Since both vectors point in the -z direction, their magnitudes add up.")
print("|B_total| = |B1| + |B2|")
print(f"|B_total| = [ mu_0 * I / (2 * pi * {r1:.0f}) ] + [ mu_0 * I / (2 * pi * {r2:.0f}) ]")
print("|B_total| = 2 * (mu_0 * I / (2 * pi))")
print("|B_total| = mu_0 * I / pi")

# Final coefficient calculation
coeff = 1 / np.pi
print(f"\nThe final result for the magnitude is (mu_0 * I) / pi, which is approximately {coeff:.4f} * mu_0 * I.")