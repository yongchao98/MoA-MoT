import math

def calculate_and_print_areas():
    """
    Calculates and prints the area for various shapes with a fixed perimeter of 1,
    demonstrating the principle of the isoperimetric problem.
    """
    # The total length of the cut is fixed at 1 meter.
    # This corresponds to the perimeter of the shape's convex hull.
    P = 1.0

    print(f"Comparing the area of different shapes with a fixed perimeter of {P} meter:")
    print("-" * 75)

    # --- A. Equilateral Triangle (Symmetry Group D3) ---
    n_triangle = 3
    # The area of a regular n-gon with perimeter P is given by the formula:
    # A = (P^2 / (4 * n)) / tan(pi / n)
    area_triangle = (P**2 / (4 * n_triangle)) / math.tan(math.pi / n_triangle)
    print("Shape: Equilateral Triangle (Symmetry Group D3)")
    print(f"Formula: Area = (P^2 / (4 * n)) / tan(pi / n)")
    print(f"Calculation: Area = ({P:.1f}**2 / (4 * {n_triangle})) / tan(pi / {n_triangle}) = {area_triangle:.6f} sq. meters")
    print()

    # --- D. Square (Symmetry Group D4) ---
    n_square = 4
    area_square = (P**2 / (4 * n_square)) / math.tan(math.pi / n_square)
    print("Shape: Square (Symmetry Group D4)")
    print(f"Formula: Area = (P^2 / (4 * n)) / tan(pi / n)")
    print(f"Calculation: Area = ({P:.1f}**2 / (4 * {n_square})) / tan(pi / {n_square}) = {area_square:.6f} sq. meters")
    print()

    # --- E. Regular Hexagon (Symmetry Group D6) ---
    n_hexagon = 6
    area_hexagon = (P**2 / (4 * n_hexagon)) / math.tan(math.pi / n_hexagon)
    print("Shape: Regular Hexagon (Symmetry Group D6)")
    print(f"Formula: Area = (P^2 / (4 * n)) / tan(pi / n)")
    print(f"Calculation: Area = ({P:.1f}**2 / (4 * {n_hexagon})) / tan(pi / {n_hexagon}) = {area_hexagon:.6f} sq. meters")
    print()

    # --- G. Circle (Symmetry Group O(2)) ---
    # The area of a circle with circumference P is given by the formula:
    # A = P^2 / (4 * pi)
    area_circle = P**2 / (4 * math.pi)
    print("Shape: Circle (Symmetry Group O(2))")
    print(f"Formula: Area = P^2 / (4 * pi)")
    print(f"Calculation: Area = {P:.1f}**2 / (4 * {math.pi:.6f}) = {area_circle:.6f} sq. meters")
    print("-" * 75)

    print("\nConclusion:")
    print("The calculations show that for a fixed perimeter, the area increases as the")
    print("number of sides of a regular polygon increases, with the circle having the")
    print("maximum possible area. Therefore, the optimal cut structure is a circle.")

if __name__ == '__main__':
    calculate_and_print_areas()