import numpy as np

def solve_square_vertices():
    """
    Calculates the vertices of a square given one point on each of its four sides.
    """
    # Let the coordinates of the marked points be P1, P2, P3, P4
    P1 = np.array([0.3511, 0.2027])
    P2 = np.array([0.6753, 0.8303])
    P3 = np.array([-0.2845, 0.9905])
    P4 = np.array([-0.128, 0.2218])

    # We assume P1, P2, P3, P4 are on consecutive sides.
    # This implies P1 and P3 are on opposite sides, and P2 and P4 are on opposite sides.
    # Vector P1-P3
    v13 = P1 - P3
    # Vector P2-P4
    v24 = P2 - P4

    # Let the direction vector of one pair of sides be u = (cos(theta), sin(theta)).
    # The side length L = |(P1-P3) . u_perp| = |(P2-P4) . u|
    # This leads to (a*ct + b*st)^2 = (c*ct + d*st)^2 where
    # (a,b) relates to v13 and u_perp, and (c,d) relates to v24 and u.
    # Let's set it up more directly. The condition is that tan(theta) must satisfy one of two equations:
    # tan(theta) = -(v13[1] + v24[0]) / (-v13[0] + v24[1])
    # tan(theta) = -(v13[1] - v24[0]) / (-v13[0] - v24[1])
    
    # Derivation from (y1-y3)cos - (x1-x3)sin = +/-((x2-x4)cos + (y2-y4)sin)
    a = P1[1] - P3[1]  # y1 - y3
    b = P3[0] - P1[0]  # x3 - x1
    c = P2[0] - P4[0]  # x2 - x4
    d = P2[1] - P4[1]  # y2 - y4
    
    # Two solutions for tan(theta)
    tan_theta1 = -(a + c) / (b + d)
    tan_theta2 = -(a - c) / (b - d)

    # Determine which solution gives the larger square
    def get_side_length(tan_theta, a_val, b_val):
        den = np.sqrt(1 + tan_theta**2)
        ct = 1 / den
        st = tan_theta / den
        return np.abs(a_val * ct + b_val * st)

    L1 = get_side_length(tan_theta1, a, b)
    L2 = get_side_length(tan_theta2, a, b)

    if L1 > L2:
        tan_theta = tan_theta1
    else:
        tan_theta = tan_theta2
        
    # Calculate sin(theta) and cos(theta) from tan(theta)
    den = np.sqrt(1 + tan_theta**2)
    ct = 1 / den
    st = tan_theta / den

    # Normal vectors for the side lines
    n1 = np.array([-st, ct]) # Normal for sides passing through P1, P3
    n2 = np.array([ct, st])  # Normal for sides passing through P2, P4

    # Constants for the line equations (n.p = const)
    c1 = np.dot(n1, P1)
    c3 = np.dot(n1, P3)
    c2 = np.dot(n2, P2)
    c4 = np.dot(n2, P4)
    
    def solve_intersection(norm1, const1, norm2, const2):
        # Solves for the intersection of two lines: norm1.X=const1 and norm2.X=const2
        M = np.array([norm1, norm2])
        K = np.array([const1, const2])
        vertex = np.linalg.solve(M, K)
        return vertex

    # The vertices are the intersections of adjacent side-lines
    # The order of points assumes P1->P2->P3->P4 is counter-clockwise on the square
    vertex_D = solve_intersection(n1, c3, n2, c4) # side DA: lines from P3 and P4
    vertex_A = solve_intersection(n1, c1, n2, c4) # side AB: lines from P1 and P4
    vertex_B = solve_intersection(n1, c1, n2, c2) # side BC: lines from P1 and P2
    vertex_C = solve_intersection(n1, c3, n2, c2) # side CD: lines from P3 and P2

    vertices = [vertex_A, vertex_B, vertex_C, vertex_D]
    
    # Sort vertices by x-coordinate
    vertices.sort(key=lambda v: v[0])
    
    # Print the results
    for v in vertices:
        print(f"({v[0]:.2f}, {v[1]:.2f})")

solve_square_vertices()