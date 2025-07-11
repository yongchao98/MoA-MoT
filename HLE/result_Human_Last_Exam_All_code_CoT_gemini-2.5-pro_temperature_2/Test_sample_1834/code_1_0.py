import math

def solve_magnetic_field():
    """
    Calculates the magnetic field at a point due to two infinite wires
    and prints the step-by-step derivation.
    """
    # Define the point of interest P
    x, y, z = 1, -1, 0

    # The problem uses a symbolic current 'I'. The constant part of the B-field
    # formula can be represented as a symbolic string.
    B_const_part = "(mu_0 * I) / (2 * pi)"

    # --- Step 1: Analyze Wire 1 (on x-axis) ---
    print("Step 1: Calculate the magnetic field from Wire 1 (on the x-axis).")
    # Perpendicular distance 'r1' from the x-axis to point P(x,y,z) is sqrt(y^2 + z^2).
    r1 = math.sqrt(y**2 + z**2)
    print(f"The point P is ({x}, {y}, {z}).")
    print(f"The perpendicular distance from the x-axis is r1 = sqrt({y}^2 + {z}^2) = {r1:.0f}.")
    print(f"The magnitude of the field is B1 = (mu_0 * I) / (2 * pi * r1).")
    # Direction is determined by the Right-Hand Rule: current in +x, point is at y<0.
    # The field points in the +z direction (or +k_hat).
    print("Using the Right-Hand Rule, the direction of B1 at P is in the positive z-direction.")
    print(f"Therefore, the vector B1 = [ {B_const_part} / {r1:.0f} ] k.")
    print("-" * 40)

    # --- Step 2: Analyze Wire 2 (on y-axis) ---
    print("Step 2: Calculate the magnetic field from Wire 2 (on the y-axis).")
    # Perpendicular distance 'r2' from the y-axis to point P(x,y,z) is sqrt(x^2 + z^2).
    r2 = math.sqrt(x**2 + z**2)
    print(f"The perpendicular distance from the y-axis is r2 = sqrt({x}^2 + {z}^2) = {r2:.0f}.")
    print(f"The magnitude of the field is B2 = (mu_0 * I) / (2 * pi * r2).")
    # Direction by Right-Hand Rule: current in +y, point is at x>0.
    # The field points in the -z direction (or -k_hat).
    print("Using the Right-Hand Rule, the direction of B2 at P is in the negative z-direction.")
    print(f"Therefore, the vector B2 = -[ {B_const_part} / {r2:.0f} ] k.")
    print("-" * 40)

    # --- Step 3: Superposition Principle ---
    print("Step 3: Apply the superposition principle to find the total magnetic field.")
    print("The total field B_total is the vector sum B_total = B1 + B2.")
    print("The final equation is:")
    print(f"B_total = [ (mu_0 * I) / (2 * pi * {r1:.0f}) ] k  +  ( -[ (mu_0 * I) / (2 * pi * {r2:.0f}) ] k )")
    print("Since the magnitudes are equal and their directions are opposite, the vectors cancel each other out.")
    print("B_total = 0")
    print("-" * 40)
    
    final_magnitude = 0
    print(f"The magnitude of the total magnetic field at ({x}, {y}, {z}) is {final_magnitude}.")

solve_magnetic_field()