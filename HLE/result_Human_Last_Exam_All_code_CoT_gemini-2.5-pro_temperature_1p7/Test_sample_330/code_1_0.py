import numpy as np

def solve_square_reconstruction():
    """
    Finds the vertices of a square given one point on each side.
    """
    # Let the four given points be P1, P2, P3, P4
    P1 = np.array([0.3511, 0.2027])
    P2 = np.array([0.6753, 0.8303])
    P3 = np.array([-0.2845, 0.9905])
    P4 = np.array([-0.128, 0.2218])

    # We assume P1, P2, P3, P4 are on consecutive sides of the square.
    # We find two potential orientations for the square. We will use one.

    # Vector from P1 to P3 and P2 to P4
    d13 = P3 - P1
    d24 = P4 - P2

    norm_d13 = np.linalg.norm(d13)
    norm_d24 = np.linalg.norm(d24)

    phi13 = np.arctan2(d13[1], d13[0])
    phi24 = np.arctan2(d24[1], d24[0])
    
    beta = phi13 - phi24
    R = norm_d13 / norm_d24
    cos_beta = np.cos(beta)
    sin_beta = np.sin(beta)
    
    # There are two solutions for the orientation angle 'alpha'.
    # We choose the second one, which leads to a valid square reconstruction.
    # tan(alpha) = (-R - sin(beta)) / cos(beta)
    tan_alpha = (-R - sin_beta) / cos_beta
    alpha = np.arctan(tan_alpha)
    
    # The orientation angle of the square's sides
    theta = alpha + phi13

    # The normal vector for the first pair of parallel sides (S1, S3)
    a = np.cos(theta)
    b = np.sin(theta)
    u = np.array([a, b])
    
    # The normal vector for the second pair of parallel sides (S2, S4)
    # This is u rotated by 90 degrees.
    v = np.array([-b, a])

    # The equations of the side lines are of the form: normal_vec . (x,y) + c = 0
    # We find the constants c1, c2, c3, c4 for each line.
    c1 = -np.dot(u, P1)
    c2 = -np.dot(v, P2)
    c3 = -np.dot(u, P3)
    c4 = -np.dot(v, P4)
    
    # Lines equations:
    # L1: a*x + b*y + c1 = 0
    # L2: -b*x + a*y + c2 = 0
    # L3: a*x + b*y + c3 = 0
    # L4: -b*x + a*y + c4 = 0

    def get_intersection(n1, const1, n2, const2):
        """Solves a system of two linear equations to find the intersection point."""
        A = np.array([n1, n2])
        B = np.array([-const1, -const2])
        try:
            intersection_point = np.linalg.solve(A, B)
            return intersection_point
        except np.linalg.LinAlgError:
            return None

    # Vertices are intersections of adjacent side lines
    V2 = get_intersection(u, c1, v, c2) # Intersection of L1 and L2
    V3 = get_intersection(v, c2, u, c3) # Intersection of L2 and L3
    V4 = get_intersection(u, c3, v, c4) # Intersection of L3 and L4
    V1 = get_intersection(v, c4, u, c1) # Intersection of L4 and L1

    vertices = [V1, V2, V3, V4]
    
    # Sort vertices based on their x-coordinate
    vertices.sort(key=lambda p: p[0])

    print("The coordinates of the vertices of the square, sorted by x-coordinate, are:")
    # "Remember in the final code you still need to output each number in the final equation!"
    # We interpret this as printing each coordinate of each vertex pair.
    for v in vertices:
        print(f"({v[0]:.2f}, {v[1]:.2f})")
        
solve_square_reconstruction()