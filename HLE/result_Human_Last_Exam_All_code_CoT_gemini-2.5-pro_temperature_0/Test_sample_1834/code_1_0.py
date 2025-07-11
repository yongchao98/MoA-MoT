import numpy as np

def solve_magnetic_field():
    """
    Calculates the magnetic field from two infinite wires at a given point.
    """
    # Define the point of interest
    point = np.array([1, -1, 0])
    print(f"Calculating the magnetic field at point P = {list(point)}\n")

    # --- Wire 1: along x-axis, current in +x direction ---
    print("--- For Wire 1 (current in +x direction) ---")
    l1_hat = np.array([1, 0, 0])
    # Find the closest point on the wire (on the x-axis) to point P
    point_on_wire1 = np.array([point[0], 0, 0])
    # Calculate the perpendicular distance vector and its magnitude
    d1_vec = point - point_on_wire1
    d1_mag = np.linalg.norm(d1_vec)
    
    # Calculate the direction of the magnetic field B1
    B1_dir = np.cross(l1_hat, d1_vec)
    # The magnitude of B1 is B0 / d1, where B0 = mu_0*I/(2*pi)
    # The vector B1 is (B0 / d1) * (B1_dir / |B1_dir|)
    # Since |B1_dir| = |l1_hat|*|d1_vec|*sin(90) = 1*d1*1 = d1,
    # B1_vec_coeff = (1/d1) * (B1_dir/d1) = B1_dir / d1**2
    B1_vec_coeff = B1_dir / (d1_mag**2)

    print(f"Perpendicular distance r1 = {d1_mag}")
    print(f"Direction vector of B1 is proportional to {list(B1_dir)}")
    print(f"The magnetic field vector B1 = (mu_0 * I / (2 * pi)) * {list(np.round(B1_vec_coeff, 5))}\n")

    # --- Wire 2: along y-axis, current in +y direction ---
    print("--- For Wire 2 (current in +y direction) ---")
    l2_hat = np.array([0, 1, 0])
    # Find the closest point on the wire (on the y-axis) to point P
    point_on_wire2 = np.array([0, point[1], 0])
    # Calculate the perpendicular distance vector and its magnitude
    d2_vec = point - point_on_wire2
    d2_mag = np.linalg.norm(d2_vec)

    # Calculate the direction of the magnetic field B2
    B2_dir = np.cross(l2_hat, d2_vec)
    B2_vec_coeff = B2_dir / (d2_mag**2)

    print(f"Perpendicular distance r2 = {d2_mag}")
    print(f"Direction vector of B2 is proportional to {list(B2_dir)}")
    print(f"The magnetic field vector B2 = (mu_0 * I / (2 * pi)) * {list(np.round(B2_vec_coeff, 5))}\n")

    # --- Total Magnetic Field ---
    print("--- Total Magnetic Field ---")
    B_total_coeff = B1_vec_coeff + B2_vec_coeff
    B_total_mag_coeff = np.linalg.norm(B_total_coeff)

    print(f"B_total = B1 + B2 = (mu_0 * I / (2 * pi)) * ({list(np.round(B1_vec_coeff, 5))} + {list(np.round(B2_vec_coeff, 5))})")
    print(f"B_total = (mu_0 * I / (2 * pi)) * {list(np.round(B_total_coeff, 5))}")
    
    # Final magnitude
    # B_total_mag = B_total_mag_coeff * (mu_0 * I / (2*pi))
    # B_total_mag = 2.0 * (mu_0 * I / (2*pi)) = mu_0 * I / pi
    print("\n--- Final Magnitude ---")
    print(f"Magnitude |B_total| = {B_total_mag_coeff} * (mu_0 * I / (2 * pi))")
    print(f"Magnitude |B_total| = {B_total_mag_coeff/2.0} * (mu_0 * I / pi)")
    print(f"So, the final magnitude is (mu_0 * I) / pi")

solve_magnetic_field()