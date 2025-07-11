import numpy as np

def solve_square_vertices():
    """
    Finds the vertices of a square given one point on each of its four sides.
    """
    # The four given points
    points = {
        'A': np.array([0.3511, 0.2027]),
        'B': np.array([0.6753, 0.8303]),
        'C': np.array([-0.2845, 0.9905]),
        'D': np.array([-0.128, 0.2218])
    }
    
    point_names = ['A', 'B', 'C', 'D']

    # There are 3 ways to pair the points on opposite sides
    pairings = [
        (('A', 'C'), ('B', 'D')),  # AC opposite, BD opposite
        (('A', 'B'), ('C', 'D')),  # AB opposite, CD opposite
        (('A', 'D'), ('B', 'C'))   # AD opposite, BC opposite
    ]

    valid_solutions = []

    for (p1_name, p3_name), (p2_name, p4_name) in pairings:
        p1 = points[p1_name]
        p3 = points[p3_name]
        p2 = points[p2_name]
        p4 = points[p4_name]
        
        # P and Q are vectors connecting opposite points
        P = p1 - p3
        Q = p2 - p4

        # The condition |u.P| = |v.Q| leads to two equations for tan(theta)
        # where u = (cos(theta), sin(theta)) and v = (-sin(theta), cos(theta))
        # 1. (P.x - Q.y)cos(t) + (P.y + Q.x)sin(t) = 0
        # 2. (P.x + Q.y)cos(t) + (P.y - Q.x)sin(t) = 0
        
        den1 = (P[1] + Q[0])
        num1 = -(P[0] - Q[1])
        
        den2 = (P[1] - Q[0])
        num2 = -(P[0] + Q[1])

        tan_thetas = []
        if den1 != 0:
            tan_thetas.append(num1 / den1)
        # The second case can also be derived by rotating the first by 90 degrees
        if den2 != 0:
            tan_thetas.append(num2 / den2)

        for tan_theta in tan_thetas:
            # Find the orientation vectors u and v
            theta = np.arctan(tan_theta)
            c, s = np.cos(theta), np.sin(theta)
            u = np.array([c, s])
            v = np.array([-s, c])

            # Determine the lines of the square for a cyclic order of points
            # Test cyclic order (p1, p2, p3, p4)
            c1 = np.dot(u, p1)
            d1 = np.dot(v, p2)
            c2 = np.dot(u, p3)
            d2 = np.dot(v, p4)

            # Inverse matrix to find intersection of lines u.x=c_val, v.x=d_val
            # x = c*c_val - s*d_val
            # y = s*c_val + c*d_val
            
            # Calculate the 4 vertices
            v1 = np.array([c*c2 - s*d2, s*c2 + c*d2]) # intersection of u.x=c2, v.x=d2
            v2 = np.array([c*c1 - s*d2, s*c1 + c*d2]) # intersection of u.x=c1, v.x=d2
            v3 = np.array([c*c1 - s*d1, s*c1 + c*d1]) # intersection of u.x=c1, v.x=d1
            v4 = np.array([c*c2 - s*d1, s*c2 + c*d1]) # intersection of u.x=c2, v.x=d1

            vertices = [v1, v2, v3, v4]
            sides = [(v2, v3), (v3, v4), (v4, v1), (v1, v2)] # sides for p1, p2, p3, p4
            
            # The original points must lie ON the segments of the square's sides
            is_valid = True
            points_on_sides = [p1, p2, p3, p4]
            for i in range(4):
                p_on_side = points_on_sides[i]
                v_start, v_end = sides[i]
                
                # Check if the point is between the vertices of the side
                # using a parameter t. If 0 <= t <= 1, it's on the segment.
                side_vec = v_end - v_start
                point_vec = p_on_side - v_start
                
                # Avoid division by zero if side is vertical or horizontal
                if np.abs(side_vec[0]) > 1e-9:
                    t = point_vec[0] / side_vec[0]
                elif np.abs(side_vec[1]) > 1e-9:
                    t = point_vec[1] / side_vec[1]
                else: # side has zero length, invalid
                    t = -1 

                # A small tolerance for floating point comparisons
                if not ( -1e-9 <= t <= 1 + 1e-9):
                    is_valid = False
                    break
            
            if is_valid:
                valid_solutions.append(vertices)

    # There should be only one valid solution
    if len(valid_solutions) == 1:
        solution_vertices = valid_solutions[0]
        # Sort vertices by x-coordinate
        sorted_vertices = sorted(solution_vertices, key=lambda v: v[0])
        
        final_coords = []
        print("The coordinates of the vertices of the square are:")
        for v in sorted_vertices:
            # Round to 2 decimal places. Handle -0.00 case.
            vx = round(v[0], 2) if abs(round(v[0], 2)) > 1e-9 else 0.0
            vy = round(v[1], 2) if abs(round(v[1], 2)) > 1e-9 else 0.0
            coord_str = f"({vx:.2f}, {vy:.2f})"
            final_coords.append(coord_str)
            print(coord_str)
            
        # Return the final answer in the required format
        return f"<<<{', '.join(final_coords)}>>>"

    return "<<<Could not find a unique solution.>>>"

final_answer = solve_square_vertices()
print(final_answer)
