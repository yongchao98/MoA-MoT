import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at point P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # The point of interest P
    x, y, z = 1, -1, 0

    print(f"This script calculates the magnetic field at point P(x, y, z) = ({x}, {y}, {z}).")
    print("The total magnetic field is the vector sum of the fields from each wire (B_total = B1 + B2).")
    print("-" * 60)

    # --- Calculation for Wire 1 (on x-axis) ---
    print("For Wire 1 (along x-axis, current in +x direction):")
    # Perpendicular distance r1 from the x-axis to P is sqrt(y^2 + z^2)
    r1 = math.sqrt(y**2 + z**2)
    print(f"The perpendicular distance r1 = sqrt({y}^2 + {z}^2) = {r1}")
    print("The magnitude is B1 = (μ₀ * I) / (2 * π * r1).")
    print("By the right-hand rule, the direction of the magnetic field B1 at P is in the +z direction.")
    print(f"So, B1_vector = (μ₀ * I) / (2 * π * {r1:.1f}) * k̂")
    print("-" * 60)

    # --- Calculation for Wire 2 (on y-axis) ---
    print("For Wire 2 (along y-axis, current in +y direction):")
    # Perpendicular distance r2 from the y-axis to P is sqrt(x^2 + z^2)
    r2 = math.sqrt(x**2 + z**2)
    print(f"The perpendicular distance r2 = sqrt({x}^2 + {z}^2) = {r2}")
    print("The magnitude is B2 = (μ₀ * I) / (2 * π * r2).")
    print("By the right-hand rule, the direction of the magnetic field B2 at P is in the -z direction.")
    print(f"So, B2_vector = - (μ₀ * I) / (2 * π * {r2:.1f}) * k̂")
    print("-" * 60)

    # --- Total Magnetic Field Calculation ---
    print("To find the total magnetic field, we sum the vectors:")
    print("B_total = B1_vector + B2_vector")
    print(f"B_total = [ (μ₀ * I) / (2 * π * {r1:.1f}) ] * k̂   -   [ (μ₀ * I) / (2 * π * {r2:.1f}) ] * k̂")

    # The magnitudes are equal and the directions are opposite.
    # Let B_term_1 = (μ₀ * I) / (2 * π * r1) and B_term_2 = (μ₀ * I) / (2 * π * r2)
    # B_total_z = B_term_1 - B_term_2
    # Since r1 = r2 = 1, B_term_1 = B_term_2.
    final_magnitude = 0.0
    
    print(f"Since r1 = {r1:.1f} and r2 = {r2:.1f}, the magnitudes of B1 and B2 are equal.")
    print("The vectors point in opposite directions (+z and -z), so they cancel each other out.")
    print(f"Final Result: |B_total| = |B1 - B2| = 0")
    print("\nThe magnitude of the magnetic field at ({}, {}, {}) is {}.".format(x, y, z, final_magnitude))

solve_magnetic_field()