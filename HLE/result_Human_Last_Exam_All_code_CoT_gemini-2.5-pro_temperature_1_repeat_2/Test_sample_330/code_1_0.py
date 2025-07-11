import numpy as np

def solve_square_vertices():
    """
    Calculates the vertices of a square given one point on each side.
    """
    # Given points
    p1 = np.array([0.3511, 0.2027])
    p2 = np.array([0.6753, 0.8303])
    p3 = np.array([-0.2845, 0.9905])
    p4 = np.array([-0.128, 0.2218])
    points = [p1, p2, p3, p4]
    
    # The problem implies P1, P2, P3, P4 are on consecutive sides.
    # This means (P1, P3) and (P2, P4) are on opposite sides.
    # Let's verify this assumption by checking the property that the distance
    # between opposite points should be equal.
    # d(P1,P3)^2 = ||P3-P1||^2
    # d(P2,P4)^2 = ||P4-P2||^2
    
    u = p1 - p3
    v = p2 - p4
    
    # This configuration is the one where ||u|| is closest to ||v||.
    # We proceed with this assumption.
    # Let u = P1 - P3 and v = P2 - P4
    ux, uy = u[0], u[1]
    vx, vy = v[0], v[1]
    
    # There are two solutions for the angle of the square's sides.
    # The tangent of the angle alpha is given by:
    tan_a1 = (vx + uy) / (ux - vy)
    tan_a2 = (uy - vx) / (ux + vy)
    
    solutions = []
    
    for i, tan_a in enumerate([tan_a1, tan_a2]):
        alpha = np.arctan(tan_a)
        sin_a = np.sin(alpha)
        cos_a = np.cos(alpha)

        # Line equations for the sides of the square
        # L1 through p1: (x-p1x)sin(a) - (y-p1y)cos(a) = 0
        # -> x*sin(a) - y*cos(a) = p1x*sin(a) - p1y*cos(a)
        # L2 through p2: (x-p2x)cos(a) + (y-p2y)sin(a) = 0
        # -> x*cos(a) + y*sin(a) = p2x*cos(a) + p2y*sin(a)
        
        c1 = p1[0] * sin_a - p1[1] * cos_a
        c2 = p2[0] * cos_a + p2[1] * sin_a
        c3 = p3[0] * sin_a - p3[1] * cos_a
        c4 = p4[0] * cos_a + p4[1] * sin_a
        
        # Matrix for solving systems of linear equations
        A_matrix = np.array([[sin_a, -cos_a], [cos_a, sin_a]])
        
        # Vertices are intersections of these lines
        # Vertex A = L1 intersect L4
        b_A = np.array([c1, c4])
        vert_A = np.linalg.solve(A_matrix, b_A)
        
        # Vertex B = L1 intersect L2
        b_B = np.array([c1, c2])
        vert_B = np.linalg.solve(A_matrix, b_B)

        # Vertex C = L3 intersect L2
        b_C = np.array([c3, c2])
        vert_C = np.linalg.solve(A_matrix, b_C)

        # Vertex D = L3 intersect L4
        b_D = np.array([c3, c4])
        vert_D = np.linalg.solve(A_matrix, b_D)
        
        side_length_sq = np.sum((vert_A - vert_B)**2)
        solutions.append({'vertices': [vert_A, vert_B, vert_C, vert_D], 'area': side_length_sq})

    # Choose the solution with the larger area
    best_solution = max(solutions, key=lambda s: s['area'])
    vertices = best_solution['vertices']
    
    # Sort vertices by x-coordinate
    vertices.sort(key=lambda v: v[0])
    
    print("The coordinates of the vertices of the square are:")
    for v in vertices:
        print(f"({v[0]:.2f},{v[1]:.2f})")

solve_square_vertices()