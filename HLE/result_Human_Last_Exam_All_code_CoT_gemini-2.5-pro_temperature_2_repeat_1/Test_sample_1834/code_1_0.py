import sympy

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Define symbols for the variables in the equations
    # I is the current, mu_0 is the permeability of free space, pi is the mathematical constant pi
    I, mu_0, pi = sympy.symbols('I mu_0 pi', real=True, positive=True)

    # Define the coordinates of the point P
    x, y, z = 1, -1, 0
    
    print(f"This script calculates the magnetic field at point P(x, y, z) = ({x}, {y}, {z}).")
    print("-" * 50)

    # --- Step 1: Calculate the magnetic field from Wire 1 (on x-axis) ---
    print("Step 1: Analyzing the field from Wire 1 (current I along +x axis)")
    
    # The perpendicular distance r1 from the wire (x-axis) to the point P is sqrt(y^2 + z^2)
    r1_val = sympy.sqrt(y**2 + z**2)
    print(f"The perpendicular distance from Wire 1 to P is r1 = sqrt({y}^2 + {z}^2) = {r1_val}")

    # Magnitude of the magnetic field B1 using B = (mu_0 * I) / (2 * pi * r)
    # The direction, determined by the right-hand rule, is in the negative z-direction.
    B1_mag = (mu_0 * I) / (2 * pi * r1_val)
    print(f"The magnitude of the magnetic field from Wire 1 is B1 = (mu_0 * I) / (2 * pi * {r1_val})")
    print("The direction of B1 at P is in the negative z-direction.")
    print("-" * 50)

    # --- Step 2: Calculate the magnetic field from Wire 2 (on y-axis) ---
    print("Step 2: Analyzing the field from Wire 2 (current I along +y axis)")

    # The perpendicular distance r2 from the wire (y-axis) to the point P is sqrt(x^2 + z^2)
    r2_val = sympy.sqrt(x**2 + z**2)
    print(f"The perpendicular distance from Wire 2 to P is r2 = sqrt({x}^2 + {z}^2) = {r2_val}")

    # Magnitude of the magnetic field B2
    # The direction, determined by the right-hand rule, is also in the negative z-direction.
    B2_mag = (mu_0 * I) / (2 * pi * r2_val)
    print(f"The magnitude of the magnetic field from Wire 2 is B2 = (mu_0 * I) / (2 * pi * {r2_val})")
    print("The direction of B2 at P is in the negative z-direction.")
    print("-" * 50)

    # --- Step 3: Calculate the total magnetic field ---
    print("Step 3: Calculating the total magnetic field using superposition")
    
    # Both B1 and B2 vectors point in the negative z-direction.
    # The total field B_total is the vector sum B1 + B2.
    # The magnitude of the total field is the sum of the individual magnitudes.
    total_magnitude_expression = B1_mag + B2_mag
    
    print(f"The total magnitude is the sum of the magnitudes of B1 and B2 because they are in the same direction.")
    print(f"|B_total| = B1 + B2 = (mu_0 * I) / (2 * pi * {r1_val}) + (mu_0 * I) / (2 * pi * {r2_val})")
    
    # Simplify the final expression
    final_magnitude = sympy.simplify(total_magnitude_expression)

    print("\nFinal simplified result:")
    print(f"The magnitude of the total magnetic field is: {final_magnitude}")

if __name__ == '__main__':
    solve_magnetic_field()