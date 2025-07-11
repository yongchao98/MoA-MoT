import math

def calculate_magnetic_field():
    """
    Calculates and explains the determination of the magnetic field at point (1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Define the point of interest
    x, y, z = 1, -1, 0

    print("--- Magnetic Field Calculation ---")
    print("This script calculates the magnetic field B at a specific point from two infinite wires.")
    print("Let the current be 'I' and the permeability of free space be 'mu_0'.")
    print(f"\nThe point of interest is P = ({x}, {y}, {z}).\n")

    # --- Wire 1 ---
    print("--- Step 1: Contribution from Wire 1 (on x-axis) ---")
    # The perpendicular distance from a point (x, y, z) to the x-axis is sqrt(y^2 + z^2)
    r1 = math.sqrt(y**2 + z**2)
    print(f"The perpendicular distance from P to Wire 1 (x-axis) is r1 = sqrt({y}^2 + {z}^2) = {r1:.1f}")
    
    print("The magnitude of the field is B1 = mu_0 * I / (2 * pi * r1).")
    print(f"B1 = mu_0 * I / (2 * pi * {r1:.1f}) = mu_0 * I / (2*pi)")
    
    print("By the right-hand rule (current in +x direction), the field at a point with a negative y-coordinate points in the +z direction.")
    print("So, the vector is B1_vec = (mu_0 * I / (2*pi)) * k_hat.\n")

    # --- Wire 2 ---
    print("--- Step 2: Contribution from Wire 2 (on y-axis) ---")
    # The perpendicular distance from a point (x, y, z) to the y-axis is sqrt(x^2 + z^2)
    r2 = math.sqrt(x**2 + z**2)
    print(f"The perpendicular distance from P to Wire 2 (y-axis) is r2 = sqrt({x}^2 + {z}^2) = {r2:.1f}")
    
    print("The magnitude of the field is B2 = mu_0 * I / (2 * pi * r2).")
    print(f"B2 = mu_0 * I / (2 * pi * {r2:.1f}) = mu_0 * I / (2*pi)")

    print("By the right-hand rule (current in +y direction), the field at a point with a positive x-coordinate points in the +z direction.")
    print("So, the vector is B2_vec = (mu_0 * I / (2*pi)) * k_hat.\n")

    # --- Total Field ---
    print("--- Step 3: Total Magnetic Field ---")
    print("The total magnetic field is the vector sum: B_total = B1_vec + B2_vec.")
    print("B_total = [ (mu_0 * I / (2*pi)) + (mu_0 * I / (2*pi)) ] * k_hat")
    print("B_total = [ 2 * (mu_0 * I / (2*pi)) ] * k_hat")
    print("B_total = (mu_0 * I / pi) * k_hat\n")

    # --- Final Magnitude ---
    print("--- Final Answer ---")
    print("The magnitude of the total magnetic field is the length of the B_total vector.")
    print("Magnitude |B| = mu_0 * I / pi")


# Run the calculation
calculate_magnetic_field()

# The final expression for the magnitude of the magnetic field.
final_answer = "mu_0 * I / pi"
# This format is requested by the prompt.
print(f"<<<{final_answer}>>>")