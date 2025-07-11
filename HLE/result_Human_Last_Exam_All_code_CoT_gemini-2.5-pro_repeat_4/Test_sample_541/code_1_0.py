import math

def solve_isoperimetric_problem():
    """
    This script demonstrates the solution to the isoperimetric problem by comparing
    the area of the convex hull for different cut shapes, each with a total length L = 1m.
    According to the Isoperimetric Theorem for Convex Hulls, the maximum area is achieved
    when the cut itself is a circle.
    """

    L = 1.0  # Total length of the cut in meters

    print(f"Analyzing shapes with a total cut length of L = {L} meter.\n")

    # --- Case G: Circle ---
    # The cut is a circle. The convex hull is the circle itself.
    # Perimeter P = 2 * pi * r = L
    r_circle = L / (2 * math.pi)
    # Area A = pi * r^2 = pi * (L / (2*pi))^2 = L^2 / (4*pi)
    area_circle = L**2 / (4 * math.pi)
    print("Case G: Cut is a Circle")
    print(f"  Equation: Area = L^2 / (4 * pi)")
    print(f"  Calculation: Area = {L**2:.1f} / (4 * {math.pi:.4f}) = {area_circle:.6f} m^2")
    print("-" * 30)

    # --- Case D: Square ---
    # The cut is a square. The convex hull is the square itself.
    # Perimeter P = 4 * s = L
    s_square = L / 4
    # Area A = s^2 = (L/4)^2 = L^2 / 16
    area_square = L**2 / 16
    print("Case D: Cut is a Square")
    print(f"  Equation: Area = L^2 / 16")
    print(f"  Calculation: Area = {L**2:.1f} / 16 = {area_square:.6f} m^2")
    print("-" * 30)
    
    # --- Case E: Regular Hexagon ---
    # The cut is a regular hexagon. The convex hull is the hexagon itself.
    # Perimeter P = 6 * s = L
    s_hexagon = L / 6
    # Area A = (3 * sqrt(3) / 2) * s^2 = (3 * sqrt(3) * L^2) / 72 = sqrt(3) * L^2 / 24
    area_hexagon = (math.sqrt(3) * L**2) / 24
    print("Case E: Cut is a Regular Hexagon")
    print(f"  Equation: Area = (sqrt(3) * L^2) / 24")
    print(f"  Calculation: Area = ({math.sqrt(3):.4f} * {L**2:.1f}) / 24 = {area_hexagon:.6f} m^2")
    print("-" * 30)

    # --- Case A: Equilateral Triangle ---
    # The cut is an equilateral triangle. The convex hull is the triangle itself.
    # Perimeter P = 3 * s = L
    s_triangle = L / 3
    # Area A = (sqrt(3) / 4) * s^2 = (sqrt(3) * L^2) / 36
    area_triangle = (math.sqrt(3) * L**2) / 36
    print("Case A: Cut is an Equilateral Triangle")
    print(f"  Equation: Area = (sqrt(3) * L^2) / 36")
    print(f"  Calculation: Area = ({math.sqrt(3):.4f} * {L**2:.1f}) / 36 = {area_triangle:.6f} m^2")
    print("-" * 30)
    
    print("\nConclusion:")
    print("The calculations show that for a fixed length, the circle encloses the maximum area.")
    print("Therefore, the optimal cut structure is a circle.")
    print("A circle has continuous rotational and reflectional symmetry, which is the symmetry group of a circle.")

solve_isoperimetric_problem()