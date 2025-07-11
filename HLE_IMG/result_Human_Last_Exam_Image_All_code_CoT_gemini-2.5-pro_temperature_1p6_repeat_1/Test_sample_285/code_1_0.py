import numpy as np

def is_in_circumcircle(p1, p2, p3, p_test, triangle_name, point_name):
    """
    Checks if p_test is inside the circumcircle of triangle (p1, p2, p3)
    using a determinant-based geometric predicate.
    It prints the matrix and determinant used for the calculation.

    Args:
        p1, p2, p3: Tuples representing the coordinates of the triangle vertices.
        p_test: A tuple representing the coordinates of the point to check.
        triangle_name: A string name for the triangle for printing.
        point_name: A string name for the test point for printing.

    Returns:
        True if p_test is inside the circumcircle, False otherwise.
    """
    ax, ay = p1
    bx, by = p2
    cx, cy = p3
    px, py = p_test

    # To ensure the test `det > 0` corresponds to "inside", the triangle
    # vertices must be ordered counter-clockwise (CCW). We check the orientation
    # and swap two vertices if the order is clockwise.
    # The orientation is determined by the sign of the 2D cross-product:
    # (bx-ax)*(cy-ay) - (by-ay)*(cx-ax)
    orientation = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    if orientation < 0:  # If clockwise, swap p2 and p3 to make it CCW.
        p2, p3 = p3, p2
        bx, by = p2
        cx, cy = p3

    # The in-circle test is based on the sign of this determinant:
    # | ax-px  ay-py  (ax-px)²+(ay-py)² |
    # | bx-px  by-py  (bx-px)²+(by-py)² |
    # | cx-px  cy-py  (cx-px)²+(cy-py)² |
    
    # We calculate each term of the matrix for clarity.
    m11 = ax - px
    m12 = ay - py
    m13 = m11**2 + m12**2
    
    m21 = bx - px
    m22 = by - py
    m23 = m21**2 + m22**2

    m31 = cx - px
    m32 = cy - py
    m33 = m31**2 + m32**2
    
    matrix = np.array([
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33]
    ])

    print(f"Checking if point {point_name}{p_test} is in the circumcircle of triangle {triangle_name}{p1, p2, p3}:")
    print("The test matrix is constructed using the coordinates of the triangle vertices and the test point:")
    print("Matrix = [[p1_x-p_x, p1_y-p_y, (p1_x-p_x)²+(p1_y-p_y)²], ...]")
    print("\nFor this test, the matrix is:")
    print(f"| {m11:8.2f} {m12:8.2f} {m13:8.2f} |")
    print(f"| {m21:8.2f} {m22:8.2f} {m23:8.2f} |")
    print(f"| {m31:8.2f} {m32:8.2f} {m33:8.2f} |")

    det_val = np.linalg.det(matrix)
    print(f"\nThe determinant of the matrix is: {det_val:.2f}")

    if det_val > 0:
        print(f"Result: The determinant is positive. Point {point_name} is INSIDE the circumcircle.")
        return True
    else:
        print(f"Result: The determinant is not positive. Point {point_name} is NOT inside the circumcircle.")
        return False

# --- Main Program ---
# 1. Define approximate coordinates for the 7 points shown in the figures.
#    (P1-P5 are outer points from top, clockwise; I1, I2 are inner points)
points_map = {
    'P1': (0, 10),    # Top vertex
    'P2': (9.5, 4),   # Top-right vertex
    'P3': (6, -8),    # Bottom-right vertex
    'P4': (-6, -8),   # Bottom-left vertex
    'P5': (-9.5, 4),  # Top-left vertex
    'I1': (-3, 1),    # Inner left point
    'I2': (3, 1)      # Inner right point
}

# 2. Focus on Figure B. The suspect triangle is T(P1, P3, P4).
triangle_vertices_keys = ['P1', 'P3', 'P4']
triangle_vertices_coords = [points_map[key] for key in triangle_vertices_keys]

# 3. The points to test are the inner points I1 and I2.
test_points_keys = ['I1', 'I2']

# 4. Perform the check for each test point.
violation_found = False
for test_key in test_points_keys:
    test_coord = points_map[test_key]
    if is_in_circumcircle(triangle_vertices_coords[0], 
                          triangle_vertices_coords[1], 
                          triangle_vertices_coords[2], 
                          test_coord, 
                          f"T({','.join(triangle_vertices_keys)})", 
                          test_key):
        violation_found = True
        print(f"\nConclusion: Triangulation B violates the empty circle property because the circumcircle of triangle T({','.join(triangle_vertices_keys)}) contains point {test_key}.")
        break # A single violation is sufficient.
    print("-" * 30)

if violation_found:
    print("\nThe analysis shows that Figure B is not a Delaunay triangulation.")
else:
    print("\nThe analysis did not find a violation for the selected triangle in Figure B.")
