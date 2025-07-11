import numpy as np

# This script calculates the magnitude of the magnetic field at a point
# due to two perpendicular infinite wires.

# --- Step 1: Define the problem setup ---
# Point of interest P
P = np.array([1, -1, 0])

# The formula for the magnetic field from a long wire is B = (mu_0 * I) / (2 * pi * r).
# We can rewrite this as B = (1 / (2 * r)) * (mu_0 * I / pi).
# The script will calculate the final coefficient for the expression (mu_0 * I / pi).

print("--- Calculating Magnetic Field at point P(1, -1, 0) ---")
print(f"Point of interest P = {list(P)}\n")

# --- Step 2: Calculate the magnetic field from Wire 1 ---
# Wire 1 lies along the x-axis with current in the +x direction.
print("--- Analyzing Wire 1 (on x-axis, current in +x direction) ---")
I1_dir = np.array([1, 0, 0])
# The perpendicular distance 'r1' is the distance from the x-axis to P.
r1 = np.sqrt(P[1]**2 + P[2]**2)
print(f"Perpendicular distance from Wire 1 to P, r1 = {r1}")

# The direction of the magnetic field B1 is found using the cross product
# of the current direction and the vector from the wire to the point.
r1_vec_comp = np.array([0, P[1], P[2]])
B1_dir_unnormalized = np.cross(I1_dir, r1_vec_comp)
B1_dir = B1_dir_unnormalized / np.linalg.norm(B1_dir_unnormalized)
print(f"Direction vector of the magnetic field from Wire 1, B1_dir = {list(B1_dir)}")

# The magnetic field vector from Wire 1, B1, in units of (mu_0 * I / pi):
coeff1 = 1 / (2 * r1)
B1_vec_coeff = coeff1 * B1_dir
print(f"Magnetic field vector from Wire 1 is B1 = {list(B1_vec_coeff)} * (mu_0 * I / pi)\n")

# --- Step 3: Calculate the magnetic field from Wire 2 ---
# Wire 2 lies along the y-axis with current in the +y direction.
print("--- Analyzing Wire 2 (on y-axis, current in +y direction) ---")
I2_dir = np.array([0, 1, 0])
# The perpendicular distance 'r2' is the distance from the y-axis to P.
r2 = np.sqrt(P[0]**2 + P[2]**2)
print(f"Perpendicular distance from Wire 2 to P, r2 = {r2}")

# The direction of the magnetic field B2 is given by the cross product.
r2_vec_comp = np.array([P[0], 0, P[2]])
B2_dir_unnormalized = np.cross(I2_dir, r2_vec_comp)
B2_dir = B2_dir_unnormalized / np.linalg.norm(B2_dir_unnormalized)
print(f"Direction vector of the magnetic field from Wire 2, B2_dir = {list(B2_dir)}")

# The magnetic field vector from Wire 2, B2, in units of (mu_0 * I / pi):
coeff2 = 1 / (2 * r2)
B2_vec_coeff = coeff2 * B2_dir
print(f"Magnetic field vector from Wire 2 is B2 = {list(B2_vec_coeff)} * (mu_0 * I / pi)\n")

# --- Step 4: Calculate the total magnetic field and its magnitude ---
print("--- Calculating Total Magnetic Field (Superposition) ---")
# The total magnetic field is the vector sum of B1 and B2.
B_total_vec_coeff = B1_vec_coeff + B2_vec_coeff
print(f"Total magnetic field vector B_total = B1 + B2 = {list(B_total_vec_coeff)} * (mu_0 * I / pi)")

# The magnitude of the total magnetic field is the norm of the total vector.
B_total_mag_coeff = np.linalg.norm(B_total_vec_coeff)
print(f"\nThe magnitude of the total field is the norm of the vector coefficient.")
print(f"|B_total| = {B_total_mag_coeff} * (mu_0 * I / pi)")

# --- Step 5: Final Answer ---
# The final equation for the magnitude is |B| = C * (mu_0 * I) / pi
print("\n--- Final Equation ---")
print("The final equation for the magnitude of the magnetic field is:")
# The numbers in the final equation as per the prompt.
print(f"|B| = ({B_total_mag_coeff}) * (mu_0 * I) / pi")
