import math

def calculate_and_compare_areas():
    """
    Calculates the area enclosed by different shapes with a fixed perimeter
    to find the one that maximizes the area.
    """
    L = 1.0  # Total length of the cut is 1 meter (perimeter/circumference)
    print("--- Comparing Areas for a Fixed Perimeter of 1.0 meter ---\n")

    # A. Equilateral Triangle
    s_triangle = L / 3
    area_triangle = (math.sqrt(3) / 4) * s_triangle**2
    print(f"If the cut forms an equilateral triangle:")
    print(f"Perimeter = 3 * Side -> Side = {L:.1f} / 3 = {s_triangle:.4f}")
    print(f"Area = (sqrt(3)/4) * Side^2 = {area_triangle:.6f}\n")

    # D. Square
    s_square = L / 4
    area_square = s_square**2
    print(f"If the cut forms a square:")
    print(f"Perimeter = 4 * Side -> Side = {L:.1f} / 4 = {s_square:.4f}")
    print(f"Area = Side^2 = {area_square:.6f}\n")

    # E. Regular Hexagon
    s_hexagon = L / 6
    area_hexagon = (3 * math.sqrt(3) / 2) * s_hexagon**2
    print(f"If the cut forms a regular hexagon:")
    print(f"Perimeter = 6 * Side -> Side = {L:.1f} / 6 = {s_hexagon:.4f}")
    print(f"Area = (3*sqrt(3)/2) * Side^2 = {area_hexagon:.6f}\n")

    # G. Circle
    r_circle = L / (2 * math.pi)
    area_circle = math.pi * r_circle**2
    print(f"If the cut forms a circle:")
    print(f"Circumference = 2*pi*Radius -> Radius = {L:.1f} / (2*pi) = {r_circle:.4f}")
    print(f"Area = pi * Radius^2 = {area_circle:.6f}\n")

    print("--- Conclusion ---")
    print(f"The area increases as the number of sides in the regular polygon increases:")
    print(f"Triangle < Square < Hexagon < ... < Circle")
    print("The maximum possible area is enclosed by the circle.")
    print("Therefore, the optimal cut structure has the symmetry of a circle.")

calculate_and_compare_areas()