import math

def get_circumcircle(p1, p2, p3):
    """
    Calculates the circumcenter and the square of the circumradius of a triangle.
    Returns (center_x, center_y, radius_sq) or None if points are collinear.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    # Using the formula for the circumcenter derived from perpendicular bisectors
    D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(D) < 1e-9:
        return None

    ux = ((x1**2 + y1**2) * (y2 - y3) + (x2**2 + y2**2) * (y3 - y1) + (x3**2 + y3**2) * (y1 - y2)) / D
    uy = ((x1**2 + y1**2) * (x3 - x2) + (x2**2 + y2**2) * (x1 - x3) + (x3**2 + y3**2) * (x2 - x1)) / D
    radius_sq = (x1 - ux)**2 + (y1 - uy)**2

    return (ux, uy, radius_sq)

def main():
    """
    Analyzes triangulation D to check for empty circle property violations.
    """
    # Estimated coordinates of the 8 points based on the image.
    points = {
        "P0_bottom_left": (0, 0),
        "P1_bottom_right": (10, 0),
        "P2_top_right": (8, 8),
        "P3_top": (5, 10),
        "P4_top_left": (2, 8),
        "P5_inner_bottom": (5, 2),
        "P6_inner_right": (6, 5),
        "P7_inner_left": (4, 5)
    }
    
    point_list = list(points.values())
    point_names = list(points.keys())

    # In triangulation D, there is a large triangle connecting non-adjacent
    # vertices of the convex hull. Let's test it.
    violating_triangle_names = ("P0_bottom_left", "P1_bottom_right", "P3_top")
    
    p1_name, p2_name, p3_name = violating_triangle_names
    p1 = points[p1_name]
    p2 = points[p2_name]
    p3 = points[p3_name]

    print(f"Analyzing Triangulation D")
    print(f"Testing triangle formed by vertices {p1_name}{p1}, {p2_name}{p2}, and {p3_name}{p3}.")

    # Calculate the circumcircle for this triangle
    circle = get_circumcircle(p1, p2, p3)
    if circle is None:
        print("The chosen vertices are collinear and do not form a triangle.")
        return

    center_x, center_y, radius_sq = circle
    print(f"\nCircumcircle Calculation:")
    print(f"Center = ({center_x:.2f}, {center_y:.2f})")
    print(f"Radius^2 = {radius_sq:.2f}")

    print("\nChecking other points against this circumcircle:")
    
    violation_found = False
    # Check all points not part of the triangle
    for i, name in enumerate(point_names):
        if name not in violating_triangle_names:
            point_to_check = point_list[i]
            
            # Calculate the squared distance from the point to the circumcenter
            dist_sq = (point_to_check[0] - center_x)**2 + (point_to_check[1] - center_y)**2
            
            is_inside = dist_sq < radius_sq
            if is_inside:
                violation_found = True

            print(f"- Point {name}{point_to_check}:")
            print(f"  Distance^2 from center = ({point_to_check[0]} - {center_x:.2f})^2 + ({point_to_check[1]} - {center_y:.2f})^2 = {dist_sq:.2f}")
            if is_inside:
                print(f"  Result: INSIDE (since {dist_sq:.2f} < {radius_sq:.2f}) -> VIOLATION")
            else:
                print(f"  Result: OUTSIDE (since {dist_sq:.2f} >= {radius_sq:.2f})")

    if violation_found:
        print("\nConclusion: At least one triangle in D violates the empty circle property.")
    else:
        print("\nConclusion: No violation found for this specific triangle.")

if __name__ == "__main__":
    main()