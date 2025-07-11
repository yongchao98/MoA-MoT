import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at point (1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Define the point P
    x, y, z = 1, -1, 0

    print("This script calculates the magnetic field at point P(1, -1, 0) from two wires.")
    print("Wire 1 is on the x-axis with current in the +x direction.")
    print("Wire 2 is on the y-axis with current in the +y direction.")
    print("The total field is the vector sum of the fields from each wire.")
    print("-" * 60)

    # --- Calculation for Wire 1 ---
    print("Field from Wire 1 (on x-axis):")
    # The perpendicular distance r1 from point P(x, y, z) to the x-axis is sqrt(y^2 + z^2)
    r1 = math.sqrt(y**2 + z**2)
    print(f"The perpendicular distance from P({x}, {y}, {z}) to the x-axis is r1 = sqrt({y}^2 + {z}^2) = {r1:.1f}")
    print("By the right-hand rule, the current in the +x direction creates a field in the -z direction at P.")
    print("So, the field vector B1 has a magnitude of (mu_0 * I) / (2 * pi * 1.0) and points in the -z direction.")
    print("-" * 60)

    # --- Calculation for Wire 2 ---
    print("Field from Wire 2 (on y-axis):")
    # The perpendicular distance r2 from point P(x, y, z) to the y-axis is sqrt(x^2 + z^2)
    r2 = math.sqrt(x**2 + z**2)
    print(f"The perpendicular distance from P({x}, {y}, {z}) to the y-axis is r2 = sqrt({x}^2 + {z}^2) = {r2:.1f}")
    print("By the right-hand rule, the current in the +y direction creates a field in the +z direction at P.")
    print("So, the field vector B2 has a magnitude of (mu_0 * I) / (2 * pi * 1.0) and points in the +z direction.")
    print("-" * 60)

    # --- Total Field Calculation ---
    print("Total Magnetic Field:")
    print("The magnitudes of B1 and B2 are equal because their distances (r1 and r2) are both 1.0.")
    print("The directions of B1 (-z) and B2 (+z) are exactly opposite.")
    print("Therefore, the two magnetic field vectors cancel each other out completely.")

    B_total_x = 0
    B_total_y = 0
    # Let B_unit = (mu_0 * I) / (2 * pi). Then B1_z = -B_unit/r1 and B2_z = B_unit/r2.
    # Since r1 = r2 = 1, the total z-component is B_unit - B_unit = 0.
    B_total_z = 0
    
    print(f"\nThe total magnetic field vector is B_total = ({B_total_x}, {B_total_y}, {B_total_z}).")

    # Calculate the final magnitude
    magnitude = math.sqrt(B_total_x**2 + B_total_y**2 + B_total_z**2)
    
    print(f"\nThe final equation for the magnitude is: |B_total| = sqrt({B_total_x}^2 + {B_total_y}^2 + {B_total_z}^2)")
    print(f"The magnitude of the magnetic field is: {magnitude}")

solve_magnetic_field()