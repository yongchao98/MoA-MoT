import numpy as np

def solve_square_vertices():
    """
    Calculates the vertices of a square given four points on its sides.
    """
    # Given points
    points = {
        'P1': np.array([0.3511, 0.2027]),
        'P2': np.array([0.6753, 0.8303]),
        'P3': np.array([-0.2845, 0.9905]),
        'P4': np.array([-0.128, 0.2218])
    }

    # Step 1: Identify the correct pairing of opposite points
    # Pair 1: (P1, P3) and (P2, P4)
    vec1 = points['P3'] - points['P1']
    vec2 = points['P4'] - points['P2']
    dot1 = np.dot(vec1, vec2)

    # Pair 2: (P1, P2) and (P3, P4)
    vec3 = points['P2'] - points['P1']
    vec4 = points['P4'] - points['P3']
    dot2 = np.dot(vec3, vec4)

    # Pair 3: (P1, P4) and (P2, P3)
    vec5 = points['P4'] - points['P1']
    vec6 = points['P3'] - points['P2']
    dot3 = np.dot(vec5, vec6)

    # The pairing with the dot product closest to 0 is the correct one.
    dots = {abs(dot1): ('P1', 'P3', 'P2', 'P4'),
            abs(dot2): ('P1', 'P2', 'P3', 'P4'),
            abs(dot3): ('P1', 'P4', 'P2', 'P3')}
    
    p1_name, p3_name, p2_name, p4_name = dots[min(dots.keys())]
    
    p1 = points[p1_name]
    p3 = points[p3_name]
    p2 = points[p2_name]
    p4 = points[p4_name]

    # Step 2 & 3: Define the directions of the square's sides
    # Vector connecting the first pair of opposite points
    v13 = p3 - p1
    # Define the normal vector for two sides based on v13
    # This will be the direction for the other two sides
    dir1 = v13
    # The perpendicular direction
    dir2 = np.array([-dir1[1], dir1[0]])

    # The normals to the side lines
    n1 = dir2
    n2 = dir1

    # Step 4: Define the four lines and find their intersections
    # Line equations are n_x*x + n_y*y = c
    # L1: through p1, perpendicular to p2p4 (normal is n2)
    c1 = np.dot(n2, p1)
    # L2: through p2, perpendicular to p1p3 (normal is n1)
    c2 = np.dot(n1, p2)
    # L3: through p3, perpendicular to p2p4 (normal is n2)
    c3 = np.dot(n2, p3)
    # L4: through p4, perpendicular to p1p3 (normal is n1)
    c4 = np.dot(n1, p4)

    # The lines are:
    # L1: n2[0]x + n2[1]y = c1
    # L2: n1[0]x + n1[1]y = c2
    # L3: n2[0]x + n2[1]y = c3
    # L4: n1[0]x + n1[1]y = c4

    # Matrix for solving systems of linear equations
    A = np.array([n2, n1])
    
    # V1 is intersection of L1 and L2
    b1 = np.array([c1, c2])
    v1 = np.linalg.solve(A, b1)

    # V2 is intersection of L2 and L3
    b2 = np.array([c3, c2])
    v2 = np.linalg.solve(A, b2)

    # V3 is intersection of L3 and L4
    b3 = np.array([c3, c4])
    v3 = np.linalg.solve(A, b3)

    # V4 is intersection of L4 and L1
    b4 = np.array([c1, c4])
    v4 = np.linalg.solve(A, b4)

    vertices = [v1, v2, v3, v4]
    
    # Step 5: Sort vertices by x-coordinate
    vertices.sort(key=lambda v: v[0])

    # Print the results
    print("The coordinates of the vertices of the square are:")
    for v in vertices:
        print(f"({v[0]:.2f}, {v[1]:.2f})")

solve_square_vertices()
<<<The coordinates of the vertices of the square are:
(-0.57, 0.76)
(0.07, -0.02)
(0.22, 1.40)
(0.85, 0.61)>>>