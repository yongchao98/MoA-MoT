import math

def calculate_and_print_area(n, shape_name):
    """
    Calculates and prints the area of the convex hull for a given number of radial cuts.
    The convex hull is a regular n-sided polygon.
    The total length of the cuts is 1.
    """
    r = 1 / n
    angle_rad = (2 * math.pi) / n
    angle_deg = 360 / n
    area = 0.5 * n * (r**2) * math.sin(angle_rad)

    print(f"For the symmetry of a {shape_name} (n={n}):")
    print(f"The cut structure consists of {n} arms, each of length 1/{n} m.")
    # The final equation part
    print(f"Area = (1/2) * {n} * (1/{n})^2 * sin(2*pi/{n})")
    print(f"Area = (1/2) * {n} * (1/{n*n}) * sin({angle_deg:.0f} degrees)")
    print(f"Area ≈ {area:.6f} square meters\n")
    return area

# --- Main execution ---
print("Calculating the area of the convex hull for different cut symmetries...\n")

# A. Equilateral Triangle (n=3)
area_triangle = calculate_and_print_area(3, "regular equilateral triangle")

# D. Square (n=4)
area_square = calculate_and_print_area(4, "square")

# E. Regular Hexagon (n=6)
area_hexagon = calculate_and_print_area(6, "regular hexagon")

if area_triangle > area_square and area_triangle > area_hexagon:
    print("Conclusion: The structure with the symmetry of an equilateral triangle maximizes the area among these options.")
else:
    print("Conclusion: The maximum area among these options is not the triangle.")
