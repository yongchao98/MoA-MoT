import math

def solve_paper_cut_problem():
    """
    Calculates and compares the convex hull area for different cut shapes
    with a total length of 1 meter.
    """
    L = 1.0
    print(f"Investigating the maximum convex hull area for cuts of total length L = {L:.1f} m.\n")

    results = {}

    # Shape G: Circle Symmetry
    # The cut is a closed circle. Perimeter = L.
    # L = 2 * pi * r  => r = L / (2 * pi)
    # Area = pi * r^2 = pi * (L / (2 * pi))^2 = L^2 / (4 * pi)
    area_circle = L**2 / (4 * math.pi)
    results['G'] = area_circle
    print("--- Shape G: Circle ---")
    print("Cut shape: A closed circle.")
    print(f"Equation: Area = {L:.1f}**2 / (4 * pi)")
    print(f"Resulting Area: {area_circle:.5f} m^2")
    print("-" * 30)

    # Function for star-shaped cuts with n arms
    def calculate_star_shape_area(n, name, key):
        # A star with n arms, each of length L/n.
        # The convex hull is a regular n-gon.
        # The distance from center to vertex (circumradius) is R = L/n.
        # Area of n-gon = n * (1/2) * R^2 * sin(2*pi/n)
        R = L / n
        angle = 2 * math.pi / n
        area = n * 0.5 * (R**2) * math.sin(angle)
        results[key] = area
        print(f"--- Shape {key}: {name} Symmetry ---")
        print(f"Cut shape: A star with {n} arms, each of length {R:.3f} m.")
        print(f"Equation: Area = {n} * 0.5 * ({L:.1f}/{n})**2 * sin(2*pi/{n})")
        print(f"Resulting Area: {area:.5f} m^2")
        print("-" * 30)

    # Shape A: Equilateral Triangle Symmetry
    calculate_star_shape_area(3, "Equilateral Triangle", "A")

    # Shape C: Isosceles Triangle Symmetry (non-equilateral)
    # Model: a 'V' shape cut. Two arms of length L/2 at 90 degrees.
    # Convex hull is a triangle with vertices (0,0), (L/2, 0), (0, L/2).
    # Area = 1/2 * base * height = 1/2 * (L/2) * (L/2) = L^2 / 8
    arm_length_v = L / 2
    area_v_shape = L**2 / 8
    results['C'] = area_v_shape
    print("--- Shape C: Isosceles Triangle Symmetry ---")
    print("Cut shape: A 'V' with two arms of length 0.5 m at 90 degrees.")
    print(f"Equation: Area = {L:.1f}**2 / 8")
    print(f"Resulting Area: {area_v_shape:.5f} m^2")
    print("-" * 30)

    # Shape D: Square Symmetry
    calculate_star_shape_area(4, "Square", "D")

    # Shape E: Regular Hexagon Symmetry
    calculate_star_shape_area(6, "Regular Hexagon", "E")

    print("\nSummary of Areas:")
    # Sorting results by area in descending order
    sorted_results = sorted(results.items(), key=lambda item: item[1], reverse=True)
    for key, area in sorted_results:
        print(f"  Shape {key}: {area:.5f} m^2")
    
    best_shape = sorted_results[0]
    print(f"\nThe best result among the tested shapes is from the one with the symmetry of an Equilateral Triangle (Shape {best_shape[0]}).")

solve_paper_cut_problem()
<<<A>>>