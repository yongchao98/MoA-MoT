import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point P(1, -1, 0)
    due to two perpendicular infinite wires.
    """
    
    # Point of interest
    P = (1, -1, 0)
    
    print("This script calculates the magnitude of the magnetic field at P(1, -1, 0).")
    print("The total field is the vector sum of the fields from each wire (Superposition Principle).\n")
    
    # --- Calculation for Wire 1 (on x-axis) ---
    print("--- Analyzing Wire 1 (on x-axis, current in +x direction) ---")
    # The closest point on the x-axis to P(1, -1, 0) is (1, 0, 0).
    r1 = 1.0
    print(f"The perpendicular distance from Wire 1 to point P is r1 = {r1}.")
    
    # Magnitude and Direction of B1
    print(f"The magnitude of the field from Wire 1 is B1 = (mu_0 * I) / (2 * pi * r1).")
    print("By the right-hand rule, the direction of the magnetic field B1 at P is in the -z direction.")
    print(f"So, B1_vector = - (mu_0 * I) / (2 * pi * {r1}) * k_hat\n")
    
    # --- Calculation for Wire 2 (on y-axis) ---
    print("--- Analyzing Wire 2 (on y-axis, current in +y direction) ---")
    # The closest point on the y-axis to P(1, -1, 0) is (0, -1, 0).
    r2 = 1.0
    print(f"The perpendicular distance from Wire 2 to point P is r2 = {r2}.")
    
    # Magnitude and Direction of B2
    print(f"The magnitude of the field from Wire 2 is B2 = (mu_0 * I) / (2 * pi * r2).")
    print("By the right-hand rule, the direction of the magnetic field B2 at P is also in the -z direction.")
    print(f"So, B2_vector = - (mu_0 * I) / (2 * pi * {r2}) * k_hat\n")
    
    # --- Total Magnetic Field ---
    print("--- Calculating the Total Magnetic Field ---")
    print("Since both B1 and B2 vectors point in the same direction (-z), we add their magnitudes.")
    print("Magnitude |B_total| = B1 + B2")
    
    print("\nSubstituting the expressions and values:")
    print(f"|B_total| = [ (mu_0 * I) / (2 * pi * {r1}) ] + [ (mu_0 * I) / (2 * pi * {r2}) ]")
    
    print("\nSimplifying the equation:")
    print(f"|B_total| = 2 * (mu_0 * I) / (2 * pi * {r1})")
    print(f"|B_total| = (mu_0 * I) / (pi * {r1})")
    
    print("\n--- Final Answer ---")
    print("The final expression for the magnitude of the magnetic field is:")
    print(f"B_magnitude = (1 * mu_0 * I) / ({r1} * pi)")
    print("Which simplifies to:")
    print("B_magnitude = (mu_0 * I) / pi")

solve_magnetic_field()