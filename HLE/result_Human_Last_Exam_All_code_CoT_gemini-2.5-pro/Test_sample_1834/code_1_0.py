import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point due to two infinite wires.
    """
    # Define the point of interest
    x, y, z = 1, -1, 0

    print(f"Calculating the magnetic field at point P(x,y,z) = ({x}, {y}, {z}).")
    print("The total field B_total is the vector sum of the fields from Wire 1 and Wire 2.")
    print("----------------------------------------------------------\n")

    # --- Analysis for Wire 1 (along x-axis, current in +x direction) ---
    print("--- Wire 1 (on x-axis) ---")
    # The perpendicular distance r1 is the distance from P to the x-axis.
    r1 = math.sqrt(y**2 + z**2)
    print(f"The perpendicular distance r1 = sqrt(y^2 + z^2) = sqrt(({y})^2 + {z}^2) = {r1}")

    print("The magnitude of the field B1 is (mu_0 * I) / (2 * pi * r1).")
    print("By the right-hand rule (current in +x, point is in -y), the field direction is -z.")
    print(f"So, B1_vector = - (mu_0 * I) / (2 * pi * {r1}) k_hat.\n")

    # --- Analysis for Wire 2 (along y-axis, current in +y direction) ---
    print("--- Wire 2 (on y-axis) ---")
    # The perpendicular distance r2 is the distance from P to the y-axis.
    r2 = math.sqrt(x**2 + z**2)
    print(f"The perpendicular distance r2 = sqrt(x^2 + z^2) = sqrt({x}^2 + {z}^2) = {r2}")

    print("The magnitude of the field B2 is (mu_0 * I) / (2 * pi * r2).")
    print("By the right-hand rule (current in +y, point is in +x), the field direction is -z.")
    print(f"So, B2_vector = - (mu_0 * I) / (2 * pi * {r2}) k_hat.\n")

    # --- Total Magnetic Field ---
    print("--- Total Magnetic Field ---")
    print("Both B1 and B2 point in the same direction (-z), so we add their magnitudes.")
    print("Magnitude |B_total| = |B1| + |B2|")
    print(f"|B_total| = (mu_0 * I) / (2 * pi * {r1}) + (mu_0 * I) / (2 * pi * {r2})")

    # Combine terms
    coeff_sum = 1/r1 + 1/r2
    print(f"|B_total| = (mu_0 * I / (2 * pi)) * (1/{r1} + 1/{r2})")
    print(f"|B_total| = (mu_0 * I / (2 * pi)) * ({coeff_sum})")
    
    # Final simplification
    final_coeff = coeff_sum / 2
    print(f"|B_total| = {final_coeff} * (mu_0 * I / pi)")

    print("\nFinal Answer:")
    print("The magnitude of the magnetic field is (mu_0 * I) / pi.")

solve_magnetic_field()