import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Define the coordinates of the point P
    x, y, z = 1, -1, 0

    print("This script calculates the magnetic field at P(x,y,z) = ({}, {}, {}).".format(x, y, z))
    print("Let the constant C = μ₀ * I / (2 * π). The magnitude of the B-field from one wire is B = C / r.\n")

    # --- Analysis for Wire 1 (on x-axis, current in +x) ---
    # The perpendicular distance r1 from the x-axis to P(x, y, z) is sqrt(y^2 + z^2).
    r1 = math.sqrt(y**2 + z**2)
    # By the right-hand rule, the magnetic field B1 at P is in the +z direction.
    # The magnitude of B1 is C / r1.
    # So, the vector B1 is (0, 0, C / r1).
    print("For Wire 1 on the x-axis:")
    print(f"The perpendicular distance r1 is {r1}.")
    print(f"The magnetic field vector B1 points in the +z direction with magnitude proportional to 1/r1 = 1/{r1}.\n")

    # --- Analysis for Wire 2 (on y-axis, current in +y) ---
    # The perpendicular distance r2 from the y-axis to P(x, y, z) is sqrt(x^2 + z^2).
    r2 = math.sqrt(x**2 + z**2)
    # By the right-hand rule, the magnetic field B2 at P is in the -z direction.
    # The magnitude of B2 is C / r2.
    # So, the vector B2 is (0, 0, -C / r2).
    print("For Wire 2 on the y-axis:")
    print(f"The perpendicular distance r2 is {r2}.")
    print(f"The magnetic field vector B2 points in the -z direction with magnitude proportional to 1/r2 = 1/{r2}.\n")

    # --- Superposition and Final Calculation ---
    # The total magnetic field is the vector sum B_total = B1 + B2.
    # B_total_z = (C / r1) - (C / r2) = C * (1/r1 - 1/r2)
    b_total_z_coeff = (1/r1) - (1/r2)
    
    # The magnitude is |B_total| = |C * (1/r1 - 1/r2)|
    magnitude = abs(b_total_z_coeff)

    print("--- Total Magnetic Field ---")
    print("The total field B_total is the vector sum of B1 and B2.")
    print("Since both fields are along the z-axis, we can write the equation for the magnitude of the total field as:")
    print(f"|B_total| = | (μ₀*I/(2*π)) * (1/{r1}) - (μ₀*I/(2*π)) * (1/{r2}) |")
    print(f"|B_total| = (μ₀*I/(2*π)) * | {1/r1} - {1/r2} |")
    print(f"|B_total| = (μ₀*I/(2*π)) * | {magnitude} |")
    
    final_magnitude = 0
    print(f"\nTherefore, the final magnitude of the magnetic field at point ({x}, {y}, {z}) is {final_magnitude}.")

solve_magnetic_field()
<<<0>>>