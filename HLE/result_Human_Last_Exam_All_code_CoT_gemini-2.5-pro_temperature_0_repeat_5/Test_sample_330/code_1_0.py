import numpy as np

def solve_square_vertices():
    """
    Calculates the vertices of a square given four points on its sides.
    """
    # Given points
    p1 = np.array([0.3511, 0.2027])
    p2 = np.array([0.6753, 0.8303])
    p3 = np.array([-0.2845, 0.9905])
    p4 = np.array([-0.128, 0.2218])

    # Step 1: Calculate vectors between opposite points
    u = p1 - p3
    v = p2 - p4

    # Step 2: Find the direction vectors for the normals of the square's sides.
    # There are two possible solutions for the square's orientation.
    # The direction vectors for the normal 'n' are derived from the condition
    # |n.u| = |n_perp.v|.
    # The two unnormalized direction vectors are:
    d1 = np.array([-u[1] - v[0], u[0] - v[1]])
    d2 = np.array([v[0] - u[1], u[0] + v[1]])

    # We choose the first solution
    normal_dir = d1
    
    # Step 3: Normalize the direction vector to get the normal vector n
    n = normal_dir / np.linalg.norm(normal_dir)
    n_perp = np.array([-n[1], n[0]])

    # Step 4: Determine the equations of the four lines of the square
    # Line equation form: ax + by = c
    c1 = np.dot(n, p1)
    c3 = np.dot(n, p3)
    c2 = np.dot(n_perp, p2)
    c4 = np.dot(n_perp, p4)

    # The lines are:
    # L1: n[0]x + n[1]y = c1
    # L3: n[0]x + n[1]y = c3
    # L2: n_perp[0]x + n_perp[1]y = c2
    # L4: n_perp[0]x + n_perp[1]y = c4

    # Step 5: Calculate the vertices by finding the intersection of these lines
    # V1 = L1 intersect L4
    # V2 = L1 intersect L2
    # V3 = L2 intersect L3
    # V4 = L3 intersect L4
    
    A14 = np.array([n, n_perp])
    b14 = np.array([c1, c4])
    v1 = np.linalg.solve(A14, b14)

    A12 = np.array([n, n_perp])
    b12 = np.array([c1, c2])
    v2 = np.linalg.solve(A12, b12)

    A23 = np.array([n_perp, n])
    b23 = np.array([c2, c3])
    v3 = np.linalg.solve(A23, b23)

    A34 = np.array([n, n_perp])
    b34 = np.array([c3, c4])
    v4 = np.linalg.solve(A34, b34)

    vertices = [v1, v2, v3, v4]

    # Step 6: Sort vertices by x-coordinate
    vertices.sort(key=lambda p: p[0])

    print("The coordinates of the vertices of the square are:")
    for vertex in vertices:
        print(f"({vertex[0]:.2f}, {vertex[1]:.2f})")

solve_square_vertices()