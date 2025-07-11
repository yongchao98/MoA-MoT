import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at point P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Point of interest P(x, y, z)
    x, y, z = 1, -1, 0

    print(f"The problem is to find the magnetic field at the point P({x}, {y}, {z}).")
    print("The total field is the vector sum of the fields from each wire (Principle of Superposition).")
    print("B_total = B_wire1 + B_wire2")
    print("-" * 50)

    # --- Wire 1: Along the x-axis, current in +x direction ---
    print("Step 1: Analyzing the field from Wire 1 (on the x-axis)")

    # The perpendicular distance from a point (x, y, z) to the x-axis is sqrt(y^2 + z^2)
    r1 = math.sqrt(y**2 + z**2)
    print(f"The perpendicular distance from P({x}, {y}, {z}) to the x-axis is r1 = sqrt({y}^2 + {z}^2) = {r1}")

    # Direction using the Right-Hand Rule:
    # Thumb in +x direction. At a point with y<0, fingers curl out of the xy-plane (+z direction).
    print("Using the Right-Hand Rule, the field B1 at P points in the positive z-direction.")

    print(f"The magnetic field vector from Wire 1 is B1 = [(μ₀ * I) / (2 * π * r1)] k̂")
    print("-" * 50)

    # --- Wire 2: Along the y-axis, current in +y direction ---
    print("Step 2: Analyzing the field from Wire 2 (on the y-axis)")

    # The perpendicular distance from a point (x, y, z) to the y-axis is sqrt(x^2 + z^2)
    r2 = math.sqrt(x**2 + z**2)
    print(f"The perpendicular distance from P({x}, {y}, {z}) to the y-axis is r2 = sqrt({x}^2 + {z}^2) = {r2}")

    # Direction using the Right-Hand Rule:
    # Thumb in +y direction. At a point with x>0, fingers curl out of the xy-plane (+z direction).
    print("Using the Right-Hand Rule, the field B2 at P also points in the positive z-direction.")
    
    print(f"The magnetic field vector from Wire 2 is B2 = [(μ₀ * I) / (2 * π * r2)] k̂")
    print("-" * 50)

    # --- Total Field ---
    print("Step 3: Calculating the total magnetic field")
    print("Since both B1 and B2 point in the same direction (+z), we add their magnitudes.")
    print("Magnitude |B_total| = |B1| + |B2|")
    
    # Display the final equation with numbers
    print("\nFinal Calculation:")
    print(f"Magnitude = (μ₀ * I) / (2 * π * {r1}) + (μ₀ * I) / (2 * π * {r2})")
    print("= (μ₀ * I) / (2 * π) + (μ₀ * I) / (2 * π)")
    print("= 2 * (μ₀ * I) / (2 * π)")
    print("= (μ₀ * I) / π")

solve_magnetic_field()