import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at point P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Define the point of interest P(x, y, z)
    P = (1, -1, 0)
    x, y, z = P
    print(f"Calculating the magnetic field at point P(x,y,z) = ({x}, {y}, {z}).\n")

    # The problem is symbolic, so let's represent the constant term (μ₀ * I) / (2 * π) as a symbol B_factor.
    # The magnetic field magnitude from a wire is B = B_factor / r.
    symbolic_B_factor = "μ₀*I/(2*π)"

    # --- Calculation for Wire 1 (on x-axis, current in +x direction) ---
    print("--- Wire 1 (on x-axis) ---")
    # The perpendicular distance r1 from point P to the x-axis is sqrt(y^2 + z^2)
    r1 = math.sqrt(y**2 + z**2)
    print(f"Perpendicular distance r1 from P to x-axis = sqrt(({y})^2 + {z}^2) = {r1}")
    # Using the right-hand rule (current in +x), the magnetic field at P(1, -1, 0) points in the -z direction.
    print("Direction of B1 field is along the negative z-axis.")
    # The magnetic field vector from wire 1 is B1 = (0, 0, -B_factor/r1)
    print(f"Symbolic B1 vector = (0, 0, -{symbolic_B_factor} / r1) = (0, 0, -{symbolic_B_factor} / {r1})\n")

    # --- Calculation for Wire 2 (on y-axis, current in +y direction) ---
    print("--- Wire 2 (on y-axis) ---")
    # The perpendicular distance r2 from point P to the y-axis is sqrt(x^2 + z^2)
    r2 = math.sqrt(x**2 + z**2)
    print(f"Perpendicular distance r2 from P to y-axis = sqrt({x}^2 + {z}^2) = {r2}")
    # Using the right-hand rule (current in +y), the magnetic field at P(1, -1, 0) points in the +z direction.
    print("Direction of B2 field is along the positive z-axis.")
    # The magnetic field vector from wire 2 is B2 = (0, 0, +B_factor/r2)
    print(f"Symbolic B2 vector = (0, 0, +{symbolic_B_factor} / r2) = (0, 0, +{symbolic_B_factor} / {r2})\n")

    # --- Total Magnetic Field (by Superposition) ---
    print("--- Total Magnetic Field ---")
    # The total magnetic field is the vector sum B_total = B1 + B2.
    print("B_total = B1 + B2")
    B_total_x = 0
    B_total_y = 0
    # The z-components from each wire are equal in magnitude and opposite in direction.
    print(f"B_total_z = (z-component of B2) + (z-component of B1)")
    print(f"B_total_z = (+{symbolic_B_factor}/{r2}) - ({symbolic_B_factor}/{r1})")
    # Since r1 = 1 and r2 = 1, the two terms cancel out.
    B_total_z = 0
    print(f"B_total_z = 0")
    print(f"Total B field vector B_total = ({B_total_x}, {B_total_y}, {B_total_z})\n")

    # --- Magnitude of the Total Magnetic Field ---
    print("--- Final Magnitude Calculation ---")
    # Magnitude = sqrt(Bx^2 + By^2 + Bz^2)
    magnitude = math.sqrt(B_total_x**2 + B_total_y**2 + B_total_z**2)
    print("The final equation for the magnitude is |B_total| = sqrt(B_total_x^2 + B_total_y^2 + B_total_z^2)")
    print(f"Magnitude |B_total| = sqrt({B_total_x}^2 + {B_total_y}^2 + {B_total_z}^2)")
    print(f"Magnitude |B_total| = sqrt({int(B_total_x**2)} + {int(B_total_y**2)} + {int(B_total_z**2)})")
    print(f"Magnitude |B_total| = {magnitude}")

# Execute the function to print the solution steps.
solve_magnetic_field()