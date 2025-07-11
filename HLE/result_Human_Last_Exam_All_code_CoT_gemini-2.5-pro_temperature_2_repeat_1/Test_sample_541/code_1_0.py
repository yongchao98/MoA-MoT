import math

def calculate_n_star_area(n):
    """
    Calculates the convex hull area for a symmetric n-star of total length 1.
    The convex hull is a regular n-gon.
    
    Args:
        n (int): The number of spokes in the star.

    Returns:
        float: The area of the convex hull.
    """
    if n < 2:
        return 0
    # For a total length of 1, each of the n spokes has length R = 1/n.
    # The area of the resulting regular n-gon is given by the formula:
    # Area = (n * R^2 / 2) * sin(2*pi/n)
    # Substituting R = 1/n, we get:
    # Area = (n * (1/n)^2 / 2) * sin(2*pi/n) = (1 / (2*n)) * sin(2*pi/n)
    area = (1 / (2 * n)) * math.sin(2 * math.pi / n)
    return area

def calculate_circular_cut_area():
    """
    Calculates the area for a circular cut of circumference 1.
    
    Returns:
        float: The area of the circle.
    """
    # For a perimeter P=1, the area of a circle is A = P^2 / (4*pi).
    area = 1 / (4 * math.pi)
    return area

def main():
    """
    Compares convex hull areas for different symmetric cut structures.
    """
    print("Investigating the area of the convex hull for different cut structures of total length 1.\n")
    
    results = {}

    # A. That of an equilateral triangle (Symmetry group D_3)
    n_triangle = 3
    area_triangle = calculate_n_star_area(n_triangle)
    results['A'] = area_triangle
    print("Symmetry: Equilateral Triangle (D3)")
    print(f"Model: A 3-spoke star with each spoke of length 1/{n_triangle}.")
    print(f"The convex hull is an equilateral triangle.")
    print(f"Area Formula: (1 / (2 * {n_triangle})) * sin(2*pi/{n_triangle})")
    print(f"Resulting Area: {area_triangle:.6f}\n")
    
    # D. That of a square (Symmetry group D_4)
    n_square = 4
    area_square = calculate_n_star_area(n_square)
    results['D'] = area_square
    print("Symmetry: Square (D4)")
    print(f"Model: A 4-spoke star with each spoke of length 1/{n_square}.")
    print(f"The convex hull is a square.")
    print(f"Area Formula: (1 / (2 * {n_square})) * sin(2*pi/{n_square})")
    print(f"Resulting Area: {area_square:.6f}\n")

    # E. That of a regular hexagon (Symmetry group D_6)
    n_hexagon = 6
    area_hexagon = calculate_n_star_area(n_hexagon)
    results['E'] = area_hexagon
    print("Symmetry: Regular Hexagon (D6)")
    print(f"Model: A 6-spoke star with each spoke of length 1/{n_hexagon}.")
    print(f"The convex hull is a regular hexagon.")
    print(f"Area Formula: (1 / (2 * {n_hexagon})) * sin(2*pi/{n_hexagon})")
    print(f"Resulting Area: {area_hexagon:.6f}\n")
    
    # G. That of a circle (Symmetry group O(2))
    area_circle = calculate_circular_cut_area()
    results['G'] = area_circle
    print("Symmetry: Circle (O(2))")
    print("Model: A single closed circular cut of length 1.")
    print("The convex hull is the circle itself.")
    print("Area Formula: 1 / (4 * pi)")
    print(f"Resulting Area: {area_circle:.6f}\n")
    
    # --- Conclusion ---
    print("-" * 40)
    print("Comparison Summary:")
    for option, area in results.items():
        print(f"  Option {option}: Area = {area:.6f}")
    
    # Find the option with the maximum area
    best_option = max(results, key=results.get)
    print(f"\nThe maximum area is achieved by the structure with the symmetry of an equilateral triangle (Option {best_option}).")
    print("-" * 40)


if __name__ == "__main__":
    main()