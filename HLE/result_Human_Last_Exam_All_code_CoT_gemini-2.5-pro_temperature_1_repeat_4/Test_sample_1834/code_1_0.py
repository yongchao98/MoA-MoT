import math

# The problem is symbolic with respect to the current I and the constant mu_0.
# The magnetic field from an infinite wire is B = (mu_0 * I) / (2 * pi * r).
# We can calculate the final magnitude as a coefficient multiplied by (mu_0 * I / (2 * pi)).

# The point of interest P
P = (1, -1, 0)
print(f"Calculating the magnetic field at point P = {P}\n")

# --- Contribution from Wire 1 (along x-axis, current in +x) ---
print("--- For Wire 1 (along x-axis) ---")
# The perpendicular distance r1 from the x-axis to a point (x,y,z) is sqrt(y^2 + z^2).
r1 = math.sqrt(P[1]**2 + P[2]**2)
# The magnitude of B1 is proportional to 1/r1.
B1_mag_coeff = 1 / r1
# The direction is found using the right-hand rule. Current is in +x.
# At P(1, -1, 0), which is in the y<0 region, the field points in the +z direction.
B1_vec_coeff = (0, 0, B1_mag_coeff)
print(f"The perpendicular distance from Wire 1 to P is r1 = sqrt(({P[1]})**2 + ({P[2]})**2) = {r1:.1f}")
print(f"The direction of B1 at P is in the +z direction.")
print(f"The vector B1 is proportional to (0, 0, {B1_vec_coeff[2]:.1f})\n")


# --- Contribution from Wire 2 (along y-axis, current in +y) ---
print("--- For Wire 2 (along y-axis) ---")
# The perpendicular distance r2 from the y-axis to a point (x,y,z) is sqrt(x^2 + z^2).
r2 = math.sqrt(P[0]**2 + P[2]**2)
# The magnitude of B2 is proportional to 1/r2.
B2_mag_coeff = 1 / r2
# The direction is found using the right-hand rule. Current is in +y.
# At P(1, -1, 0), which is in the x>0 region, the field points in the -z direction.
B2_vec_coeff = (0, 0, -B2_mag_coeff)
print(f"The perpendicular distance from Wire 2 to P is r2 = sqrt(({P[0]})**2 + ({P[2]})**2) = {r2:.1f}")
print(f"The direction of B2 at P is in the -z direction.")
print(f"The vector B2 is proportional to (0, 0, {B2_vec_coeff[2]:.1f})\n")

# --- Total Field Calculation ---
print("--- Total Magnetic Field and Magnitude ---")
# The total field is the vector sum B_total = B1 + B2.
# We add the coefficient vectors to find the total coefficient vector.
B_total_vec_coeff_x = B1_vec_coeff[0] + B2_vec_coeff[0]
B_total_vec_coeff_y = B1_vec_coeff[1] + B2_vec_coeff[1]
B_total_vec_coeff_z = B1_vec_coeff[2] + B2_vec_coeff[2]

# The final magnitude is the magnitude of the total vector.
# Let C = (mu_0 * I / (2 * pi)).
# |B_total| = C * |B1_vec_coeff + B2_vec_coeff|
print("The final magnitude equation is |B_total| = (mu_0 * I / (2 * pi)) * |(B1_vector_coeffs) + (B2_vector_coeffs)|")
print(f"|B_total| = (mu_0 * I / (2 * pi)) * |(0, 0, {B1_vec_coeff[2]:.1f}) + (0, 0, {B2_vec_coeff[2]:.1f})|")
print(f"|B_total| = (mu_0 * I / (2 * pi)) * |(0, 0, {B_total_vec_coeff_z:.1f})|")

final_magnitude_coeff = math.sqrt(B_total_vec_coeff_x**2 + B_total_vec_coeff_y**2 + B_total_vec_coeff_z**2)
print(f"\nAs the two magnetic field vectors are equal in magnitude and opposite in direction, they cancel each other out.")
print(f"The magnitude of the total magnetic field is {final_magnitude_coeff:.1f}.")
