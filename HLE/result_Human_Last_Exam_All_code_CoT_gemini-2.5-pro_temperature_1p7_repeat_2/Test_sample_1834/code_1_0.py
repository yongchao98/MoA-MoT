import numpy as np

def solve_magnetic_field():
    """
    Calculates the magnetic field at a point due to two infinite wires.
    """
    # For calculation purposes, let's represent the constant term mu_0 * I / (2 * pi)
    # as a single factor, B_factor. We can set it to 1.0 since we are interested
    # in the final magnitude, and the two individual magnitudes are identical.
    B_factor = 1.0

    # The point where we want to find the magnetic field
    P = np.array([1, -1, 0])

    print("--- Calculating Magnetic Field at P(1, -1, 0) ---")
    print(f"Let the constant term mu_0*I/(2*pi) be represented as B_unit = {B_factor:.1f}.\n")

    # --- Wire 1: along x-axis, current in +x direction ---
    r1 = np.sqrt(P[1]**2 + P[2]**2)
    B1_mag = B_factor / r1
    # By RHR, direction is +z for y < 0
    B1_vec = np.array([0, 0, B1_mag])

    print("Step 1: Calculate the magnetic field from Wire 1 (on x-axis)")
    print(f"Perpendicular distance r1 = sqrt(y^2 + z^2) = sqrt(({P[1]})**2 + {P[2]}**2) = {r1:.2f}")
    print(f"Magnitude B1 = B_unit / r1 = {B_factor:.1f} / {r1:.2f} = {B1_mag:.2f} * B_unit")
    print(f"Direction is +z. So, B1 vector = <{B1_vec[0]:.2f}, {B1_vec[1]:.2f}, {B1_vec[2]:.2f}>\n")

    # --- Wire 2: along y-axis, current in +y direction ---
    r2 = np.sqrt(P[0]**2 + P[2]**2)
    B2_mag = B_factor / r2
    # By RHR, direction is -z for x > 0
    B2_vec = np.array([0, 0, -B2_mag])

    print("Step 2: Calculate the magnetic field from Wire 2 (on y-axis)")
    print(f"Perpendicular distance r2 = sqrt(x^2 + z^2) = sqrt({P[0]}**2 + {P[2]}**2) = {r2:.2f}")
    print(f"Magnitude B2 = B_unit / r2 = {B_factor:.1f} / {r2:.2f} = {B2_mag:.2f} * B_unit")
    print(f"Direction is -z. So, B2 vector = <{B2_vec[0]:.2f}, {B2_vec[1]:.2f}, {B2_vec[2]:.2f}>\n")

    # --- Superposition Principle: Add the vectors ---
    B_total_vec = B1_vec + B2_vec
    B_total_mag = np.linalg.norm(B_total_vec)

    print("Step 3: Sum the fields using the superposition principle")
    print(f"B_total = B1 + B2 = <{B1_vec[0]:.2f}, {B1_vec[1]:.2f}, {B1_vec[2]:.2f}> + <{B2_vec[0]:.2f}, {B2_vec[1]:.2f}, {B2_vec[2]:.2f}>")
    print(f"B_total_x = {B1_vec[0]:.2f} + {B2_vec[0]:.2f} = {B_total_vec[0]:.2f}")
    print(f"B_total_y = {B1_vec[1]:.2f} + {B2_vec[1]:.2f} = {B_total_vec[1]:.2f}")
    print(f"B_total_z = {B1_vec[2]:.2f} + {B2_vec[2]:.2f} = {B_total_vec[2]:.2f}\n")


    print("Step 4: Calculate the final magnitude of the total magnetic field")
    print(f"|B_total| = sqrt( (Bx_total)^2 + (By_total)^2 + (Bz_total)^2 )")
    print(f"|B_total| = sqrt( ({B_total_vec[0]:.2f})^2 + ({B_total_vec[1]:.2f})^2 + ({B_total_vec[2]:.2f})^2 )")
    print(f"|B_total| = {B_total_mag:.2f} (in units of mu_0*I/(2*pi))")
    print("\nSince the two field vectors are equal in magnitude and opposite in direction, they cancel out.")
    print(f"The final magnitude is exactly {int(B_total_mag)}.")


solve_magnetic_field()