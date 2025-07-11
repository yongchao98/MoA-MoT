import math

def calculate_areas():
    """
    Calculates the area of different shapes with a fixed perimeter of 1 meter
    and prints the results, including the formulas used.
    """
    # Define the total length of the cut (perimeter)
    L = 1.0

    print(f"Comparing areas for shapes with a perimeter of {L} meter:\n")

    # --- Equilateral Triangle ---
    s_triangle = L / 3.0
    area_triangle = (math.sqrt(3.0) / 4.0) * (s_triangle ** 2)
    print("A. Equilateral Triangle (Symmetry Group D3)")
    print(f"   Formula: Area = (sqrt(3)/4) * (L/3)^2")
    print(f"   Calculation: Area = (sqrt(3)/4) * ({s_triangle:.4f})^2 = {area_triangle:.6f} sq. meters\n")

    # --- Square ---
    s_square = L / 4.0
    area_square = s_square ** 2
    print("D. Square (Symmetry Group D4)")
    print(f"   Formula: Area = (L/4)^2")
    print(f"   Calculation: Area = ({s_square:.4f})^2 = {area_square:.6f} sq. meters\n")

    # --- Regular Hexagon ---
    s_hexagon = L / 6.0
    area_hexagon = (3.0 * math.sqrt(3.0) / 2.0) * (s_hexagon ** 2)
    print("E. Regular Hexagon (Symmetry Group D6)")
    print(f"   Formula: Area = (3*sqrt(3)/2) * (L/6)^2")
    print(f"   Calculation: Area = (3*sqrt(3)/2) * ({s_hexagon:.4f})^2 = {area_hexagon:.6f} sq. meters\n")

    # --- Circle ---
    # Area A = pi * r^2 and Perimeter L = 2 * pi * r. So r = L / (2*pi).
    # Substituting r gives A = pi * (L / (2*pi))^2 = pi * L^2 / (4*pi^2) = L^2 / (4*pi)
    area_circle = (L ** 2) / (4.0 * math.pi)
    print("G. Circle (Symmetry Group O(2))")
    print(f"   Formula: Area = L^2 / (4*pi)")
    print(f"   Calculation: Area = {L**2:.1f} / (4*pi) = {area_circle:.6f} sq. meters\n")

    print("Conclusion: The calculations confirm that the circle encloses the maximum area.")

if __name__ == '__main__':
    calculate_areas()