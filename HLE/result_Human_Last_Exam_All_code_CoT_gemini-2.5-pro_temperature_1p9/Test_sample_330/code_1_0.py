import numpy as np

def solve_square_vertices():
    """
    Calculates the vertices of a square given one point on each side.
    """
    # Coordinates of the four marked points
    P = np.array([
        [0.3511, 0.2027],  # P1
        [0.6753, 0.8303],  # P2
        [-0.2845, 0.9905], # P3
        [-0.128, 0.2218]   # P4
    ])

    # There are 3 ways to pair the points on opposite sides.
    # We'll test the pairing (P1, P3) and (P2, P4). This assumes the points
    # P1, P2, P3, P4 lie on consecutive sides of the square.
    p_indices = [0, 1, 2, 3] # P1->P2->P3->P4 on sides AB, BC, CD, DA
    u_pair_indices = (0, 2)  # Points on one pair of opposite sides
    v_pair_indices = (1, 3)  # Points on the other pair
    
    P_ordered = P[p_indices]
    
    # Create vectors u and v from the pairs of points.
    u = P[u_pair_indices[0]] - P[u_pair_indices[1]]
    v = P[v_pair_indices[0]] - P[v_pair_indices[1]]
    
    # Rotate v by 90 degrees clockwise to use in the geometric condition.
    v_rot = np.array([v[1], -v[0]])
    
    # The condition for the square side directions leads to two solutions for the normal vector.
    # The normal vector n must be perpendicular to either (u - v_rot) or (u + v_rot).
    w1 = u - v_rot
    w2 = u + v_rot
    
    possible_normals = [
        np.array([-w1[1], w1[0]]), # Perpendicular to w1
        np.array([-w2[1], w2[0]])  # Perpendicular to w2
    ]
    
    final_vertices = None
    
    for n_unnormalized in possible_normals:
        # Avoid division by zero if a normal is a zero vector
        if np.linalg.norm(n_unnormalized) < 1e-9:
            continue
            
        # Normalize the vector to get n = (cos(theta), sin(theta))
        n = n_unnormalized / np.linalg.norm(n_unnormalized)
        # Get the perpendicular normal for the other pair of sides
        n_perp = np.array([-n[1], n[0]])

        # The lines are defined by n . (x,y) - d = 0.
        # Calculate the constant 'd' for each of the four lines.
        # Each line passes through its corresponding point.
        d = np.array([
            np.dot(n,      P_ordered[0]), # Line through P1 (side AB)
            np.dot(n_perp, P_ordered[1]), # Line through P2 (side BC)
            np.dot(n,      P_ordered[2]), # Line through P3 (side CD)
            np.dot(n_perp, P_ordered[3])  # Line through P4 (side DA)
        ])

        # The vertices are intersections of these lines.
        # The matrix for the system of equations for two perpendicular lines is [[n], [n_perp]].
        A_matrix = np.array([n, n_perp])
        
        # Vertex A is the intersection of lines DA and AB (lines for P4 and P1)
        V_A = np.linalg.solve(A_matrix, np.array([d[0], d[3]]))
        # Vertex B is the intersection of lines AB and BC (lines for P1 and P2)
        V_B = np.linalg.solve(A_matrix, np.array([d[0], d[1]]))
        # Vertex C is the intersection of lines BC and CD (lines for P2 and P3)
        V_C = np.linalg.solve(A_matrix, np.array([d[2], d[1]]))
        # Vertex D is the intersection of lines CD and DA (lines for P3 and P4)
        V_D = np.linalg.solve(A_matrix, np.array([d[2], d[3]]))
        
        vertices = np.array([V_A, V_B, V_C, V_D])
        # The sides connect the vertices in order A->B->C->D->A
        sides = [(V_A, V_B), (V_B, V_C), (V_C, V_D), (V_D, V_A)]
        
        # Check if this solution is valid: all given points must lie on the side *segments*.
        is_valid = True
        for i in range(4):
            p = P_ordered[i]
            v1, v2 = sides[i]
            # Use a simple bounding box check for segment inclusion.
            is_on_segment = (min(v1[0], v2[0]) - 1e-6 <= p[0] <= max(v1[0], v2[0]) + 1e-6 and
                             min(v1[1], v2[1]) - 1e-6 <= p[1] <= max(v1[1], v2[1]) + 1e-6)
            if not is_on_segment:
                is_valid = False
                break
        
        if is_valid:
            final_vertices = vertices
            break
            
    if final_vertices is not None:
        # Sort vertices by their x-coordinate
        sorted_vertices = sorted(final_vertices, key=lambda v: v[0])
        
        # Print the results
        for vertex in sorted_vertices:
            print(f"({vertex[0]:.2f}, {vertex[1]:.2f})")

solve_square_vertices()