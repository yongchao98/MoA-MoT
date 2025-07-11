import numpy as np

def solve_square(points):
    """
    Solves for the vertices of a square given four points, one on each side.
    The main idea is to test all 3 unique pairings of points on the square's sides.
    
    For a given consecutive ordering of points Q1, Q2, Q3, Q4 on the sides:
    Let the side direction be (c, s) = (cos(theta), sin(theta)).
    The normals to the sides are (s, -c) and (c, s).
    
    The four lines of the square can be written as:
    L1: s*x - c*y = s*x1 - c*y1
    L2: c*x + s*y = c*x2 + s*y2
    L3: s*x - c*y = s*x3 - c*y3
    L4: c*x + s*y = c*x4 + s*y4
    
    The distance between L1 and L3 is S (side length), so |(s*x3-c*y3) - (s*x1-c*y1)| = S.
    The distance between L2 and L4 is S, so |(c*x2+s*y2) - (c*x4+s*y4)| = S.
    
    This gives two equations for S, which can be solved for theta.
    tan(theta) = s/c = (x2-x4+y3-y1) / (x3-x1-y2-y4)
    """
    
    # The three distinct permutations of opposite points.
    # We list one consecutive order for each case.
    permutations = [
        [0, 1, 2, 3],  # P1-P2-P3-P4
        [0, 2, 1, 3],  # P1-P3-P2-P4
        [0, 1, 3, 2]   # P1-P2-P4-P3
    ]

    for p_indices in permutations:
        q_points = [points[i] for i in p_indices]
        q1, q2, q3, q4 = q_points

        # Case 1: Normals (s, -c) and (c, s)
        num = (q2[0] - q4[0]) + (q3[1] - q1[1])
        den = (q3[0] - q1[0]) - (q2[1] - q4[1])

        if abs(den) < 1e-9:
            continue
            
        # Two possible angles from tan, differing by 180 degrees
        for sign in [1, -1]:
            theta = np.arctan(num / den)
            c, s = np.cos(theta), np.sin(theta)
            c, s = sign * c, sign * s
            
            side_len = s * (q3[0] - q1[0]) - c * (q3[1] - q1[1])

            # Side length must be positive
            if side_len < 0:
                s, c = -s, -c
                side_len = -side_len
            
            if side_len < 1e-9: continue

            # Lines equations constants
            # L1: s*x - c*y = C1
            # L2: c*x + s*y = C2
            # L3: s*x - c*y = C3
            # L4: c*x + s*y = C4
            C1 = s * q1[0] - c * q1[1]
            C2 = c * q2[0] + s * q2[1]
            # C3 and C4 based on S
            # We need to decide on which side of L1/L2 the other lines are
            # Try one orientation
            C3 = C1 + side_len
            C4 = C2 - side_len

            # Check if this orientation matches the points
            if abs(s*q3[0] - c*q3[1] - C3) > 1e-6 or abs(c*q4[0] + s*q4[1] - C4) > 1e-6:
                # Try the other orientation
                C4 = C2 + side_len
                if abs(s*q3[0] - c*q3[1] - C3) > 1e-6 or abs(c*q4[0] + s*q4[1] - C4) > 1e-6:
                   continue

            # Vertices are intersections of these lines
            # Matrix A = [[s, -c], [c, s]]. A_inv = [[s, c], [-c, s]]
            V1 = (s * C1 + c * C4, -c * C1 + s * C4) # L1, L4
            V2 = (s * C1 + c * C2, -c * C1 + s * C2) # L1, L2
            V3 = (s * C3 + c * C2, -c * C3 + s * C2) # L3, L2
            V4 = (s * C3 + c * C4, -c * C3 + s * C4) # L3, L4

            vertices = [V1, V2, V3, V4]
            # Points for segment checking
            points_on_sides = {
                0: (V1, V2), # P_on_L1 lies on V1V2 segment
                1: (V2, V3), # P_on_L2 lies on V2V3 segment
                2: (V3, V4), # P_on_L3 lies on V3V4 segment
                3: (V4, V1)  # P_on_L4 lies on V4V1 segment
            }

            valid = True
            for i in range(4):
                p = q_points[i]
                v_start, v_end = points_on_sides[i]
                # Check for collinearity and if point is on segment
                vec_side = np.array(v_end) - np.array(v_start)
                vec_point = np.array(p) - np.array(v_start)
                cross_prod = np.cross(vec_side, vec_point)
                
                if abs(cross_prod) > 1e-6: # Not collinear
                    valid = False
                    break
                
                # Check if it is on segment using dot product
                dot_prod = np.dot(vec_point, vec_side)
                if dot_prod < 0 or dot_prod > np.dot(vec_side, vec_side):
                    valid = False
                    break
            
            if valid:
                # Sort vertices by x-coordinate
                sorted_vertices = sorted(vertices, key=lambda v: v[0])
                for vx, vy in sorted_vertices:
                    print(f"({vx:.2f},{vy:.2f})")
                return

# Given coordinates
P1 = (0.3511, 0.2027)
P2 = (0.6753, 0.8303)
P3 = (-0.2845, 0.9905)
P4 = (-0.128, 0.2218)
point_list = [P1, P2, P3, P4]

solve_square(point_list)