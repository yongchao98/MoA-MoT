import numpy as np

def check_violation_in_B():
    """
    This function checks for a Delaunay violation in triangulation B.
    It defines coordinates for the vertices of a suspicious triangle and a point
    that may lie within its circumcircle. It then performs the in-circle test.
    """
    # Estimated coordinates for the points in Triangulation B
    # T: top vertex, TR: top-right, BR: bottom-right, IR: inner-right
    t_vertex = np.array([0, 10])
    tr_vertex = np.array([8, 5])
    br_vertex = np.array([8, -8])
    ir_point = np.array([3, 1])

    print("Checking Triangulation B for Delaunay violation.")
    print("Triangle vertices to test:")
    print(f"  A (Top): {t_vertex}")
    print(f"  B (Top-Right): {tr_vertex}")
    print(f"  C (Bottom-Right): {br_vertex}")
    print(f"Point to check if it's inside the circumcircle:")
    print(f"  D (Inner-Right): {ir_point}")
    print("-" * 30)

    # The in-circle test requires points to be in counter-clockwise order.
    # The order (T, TR, BR) is clockwise. Let's use (T, BR, TR) which is CCW.
    a = t_vertex
    b = br_vertex
    c = tr_vertex
    d = ir_point
    
    # Building the matrix for the determinant calculation
    matrix = np.array([
        [a[0] - d[0], a[1] - d[1], (a[0] - d[0])**2 + (a[1] - d[1])**2],
        [b[0] - d[0], b[1] - d[1], (b[0] - d[0])**2 + (b[1] - d[1])**2],
        [c[0] - d[0], c[1] - d[1], (c[0] - d[0])**2 + (c[1] - d[1])**2]
    ])

    # Calculate the determinant
    det = np.linalg.det(matrix)

    print("In-circle test calculation:")
    print("The matrix is:")
    print(matrix)
    print(f"Determinant = {det:.2f}")

    if det > 0:
        print("\nResult: The determinant is positive.")
        print("This confirms that point D (Inner-Right) is INSIDE the circumcircle of triangle ABC (T-BR-TR).")
        print("Therefore, Triangulation B violates the empty circle property of Delaunay triangulations.")
    elif det < 0:
        print("\nResult: The determinant is negative.")
        print("This means the point is OUTSIDE the circumcircle.")
    else:
        print("\nResult: The determinant is zero.")
        print("This means the four points are co-circular.")

check_violation_in_B()