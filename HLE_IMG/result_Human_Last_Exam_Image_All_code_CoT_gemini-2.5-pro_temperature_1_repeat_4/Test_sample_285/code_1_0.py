import numpy as np

def get_points():
    """
    Defines the coordinates of the 7 points.
    The points are based on a regular pentagon centered at the origin, with a radius of 10.
    V1 is the top vertex. V2-V5 are subsequent vertices counter-clockwise.
    I1 and I2 are the interior points on the y-axis.
    """
    points = {}
    # Vertices of a regular pentagon
    # Angles for a pentagon with a vertex pointing up
    angles_deg = [90, 162, 234, 306, 18] 
    vertices = ['V1', 'V5', 'V4', 'V3', 'V2'] # Match visual layout (V1-top, V2-top-right, etc.)
    for i in range(5):
        angle = np.deg2rad(angles_deg[i])
        points[vertices[i]] = (10 * np.cos(angle), 10 * np.sin(angle))
    
    # Interior points
    points['I1'] = (0, 4)  # Upper interior point
    points['I2'] = (0, -3) # Lower interior point
    return points

def is_inside_triangle(pt, v1, v2, v3):
    """
    Checks if point pt is inside the triangle defined by v1, v2, v3.
    Uses barycentric coordinates method.
    """
    def sign(p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    d1 = sign(pt, v1, v2)
    d2 = sign(pt, v2, v3)
    d3 = sign(pt, v3, v1)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    # If all signs are the same (or zero), the point is inside or on the boundary.
    return not (has_neg and has_pos)

def is_in_circumcircle(p_test, p1, p2, p3):
    """
    Checks if p_test is inside the circumcircle of the triangle (p1, p2, p3).
    This uses the in-circle test based on a determinant, which is robust.
    Assumes p1, p2, p3 are in counter-clockwise order.
    """
    ax, ay = p1
    bx, by = p2
    cx, cy = p3
    dx, dy = p_test

    matrix = np.array([
        [ax - dx, ay - dy, (ax - dx)**2 + (ay - dy)**2],
        [bx - dx, by - dy, (bx - dx)**2 + (by - dy)**2],
        [cx - dx, cy - dy, (cx - dx)**2 + (cy - dy)**2]
    ])
    
    # A positive determinant means the point is inside the circumcircle.
    return np.linalg.det(matrix) > 0

def main():
    """
    Main function to analyze the triangulations.
    """
    points = get_points()
    violators = []

    print("Analyzing Triangulations for Delaunay Property Violations:\n")

    # --- Analysis of Triangulation B ---
    print("Analysis of B:")
    # In B, we test if the large triangle (V1, V3, V4) contains other points.
    t_B = (points['V1'], points['V3'], points['V4'])
    p_check_B_1 = points['I1']
    p_check_B_2 = points['I2']
    
    is_violator_B = False
    # A point inside a triangle is always inside its circumcircle.
    if is_inside_triangle(p_check_B_1, t_B[0], t_B[1], t_B[2]):
        print(f"Violation found in B: Point I1 {p_check_B_1} is inside triangle (V1, V3, V4).")
        is_violator_B = True
    if is_inside_triangle(p_check_B_2, t_B[0], t_B[1], t_B[2]):
        print(f"Violation found in B: Point I2 {p_check_B_2} is inside triangle (V1, V3, V4).")
        is_violator_B = True
        
    if is_violator_B:
        violators.append('B')
    else:
        print("No obvious violation found in B based on the selected triangle.")
    
    print("-" * 20)

    # --- Analysis of Triangulation D ---
    print("Analysis of D:")
    # In D, we test if the circumcircle of triangle (V2, I1, V5) contains point I2.
    # The triangle vertices must be in CCW order for the test.
    # V5 -> I1 -> V2 is CCW.
    t_D_p1 = points['V5']
    t_D_p2 = points['I1']
    t_D_p3 = points['V2']
    
    p_check_D = points['I2']

    if is_in_circumcircle(p_check_D, t_D_p1, t_D_p2, t_D_p3):
        print(f"Violation found in D: Point I2 {p_check_D} is inside the circumcircle of triangle (V5, I1, V2).")
        violators.append('D')
    else:
        print("No violation found in D based on the selected triangle.")
        
    print("-" * 20)

    # --- Analysis of Triangulation C ---
    print("Analysis of C:")
    print("Triangulation C is not a valid triangulation of the point set as it leaves areas of the convex hull untriangulated. Therefore, it cannot be a Delaunay triangulation.")
    print("-" * 20)
    
    # --- Final Conclusion ---
    if violators:
        result = ",".join(sorted(violators))
        print(f"\nConclusion: The triangulations that violate the empty circle property are {result}.")
    else:
        print("\nConclusion: No violations were found among the analyzed triangulations.")

    # The final answer in the required format
    final_answer = ",".join(sorted(violators))
    print(f"\nFinal Answer: {final_answer}")


if __name__ == "__main__":
    main()