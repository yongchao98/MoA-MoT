import numpy as np

def solve_square_vertices(points):
    """
    Finds the vertices of a square given one point on each side.

    The method involves:
    1. Trying all 3 possible pairings of points on opposite sides.
    2. For each pairing, calculating the two possible slopes (m) for the square's sides.
    3. For each slope, determining the lines and vertices of the potential square.
    4. Validating if the given points lie on the *segments* of the square's sides.
    5. Returning the vertices of the first valid square found.
    """
    
    # A small tolerance for floating point comparisons
    TOL = 1e-9

    def is_on_segment(p, v1, v2):
        """Check if point p lies on the segment between v1 and v2."""
        # Check if p is collinear with v1 and v2 (already known from line construction)
        # Check if p is within the bounding box of the segment
        in_bounds = (min(v1[0], v2[0]) - TOL <= p[0] <= max(v1[0], v2[0]) + TOL and
                     min(v1[1], v2[1]) - TOL <= p[1] <= max(v1[1], v2[1]) + TOL)
        return in_bounds

    # The 3 ways to pair 4 points
    pairings = [
        ((0, 2), (1, 3)),  # (P1, P3) opposite, (P2, P4) opposite
        ((0, 1), (2, 3)),  # (P1, P2) opposite, (P3, P4) opposite
        ((0, 3), (1, 2))   # (P1, P4) opposite, (P2, P3) opposite
    ]

    for p_indices in pairings:
        p1_idx, p3_idx = p_indices[0]
        p2_idx, p4_idx = p_indices[1]
        
        p1 = points[p1_idx]
        p2 = points[p2_idx]
        p3 = points[p3_idx]
        p4 = points[p4_idx]

        # Calculate deltas for the slope formulas
        dx13 = p1[0] - p3[0]
        dy13 = p1[1] - p3[1]
        dx24 = p2[0] - p4[0]
        dy24 = p2[1] - p4[1]
        
        # Calculate the two possible slopes from the two cases of the absolute value equation
        # Case 1: dy13 - m*dx13 = dx24 + m*dy24 => m = (dy13 - dx24) / (dx13 + dy24)
        # Case 2: dy13 - m*dx13 = -(dx24 + m*dy24) => m = (dy13 + dx24) / (dx13 - dy24)
        
        slopes = []
        if abs(dx13 + dy24) > TOL:
            slopes.append((dy13 - dx24) / (dx13 + dy24))
        if abs(dx13 - dy24) > TOL:
            slopes.append((dy13 + dx24) / (dx13 - dy24))

        for m in slopes:
            if abs(m) < TOL: continue # Avoid horizontal/vertical cases for simplicity if m is near zero
            m_perp = -1.0 / m

            # Lines equations: y = mx + c  => c = y - mx
            # L1 passes through p1 with slope m
            # L3 passes through p3 with slope m
            # L2 passes through p2 with slope m_perp
            # L4 passes through p4 with slope m_perp
            c1 = p1[1] - m * p1[0]
            c3 = p3[1] - m * p3[0]
            c2 = p2[1] - m_perp * p2[0]
            c4 = p4[1] - m_perp * p4[0]

            # Intersection of two lines y=m1x+c1 and y=m2x+c2 is x=(c2-c1)/(m1-m2)
            # Denominator for all intersections will be (m - m_perp)
            m_diff = m - m_perp
            
            # Find vertices
            # A = L1 intersect L4
            x_A = (c4 - c1) / m_diff
            y_A = m * x_A + c1
            A = (x_A, y_A)
            
            # B = L1 intersect L2
            x_B = (c2 - c1) / m_diff
            y_B = m * x_B + c1
            B = (x_B, y_B)

            # C = L3 intersect L2
            x_C = (c2 - c3) / m_diff
            y_C = m * x_C + c3
            C = (x_C, y_C)

            # D = L3 intersect L4
            x_D = (c4 - c3) / m_diff
            y_D = m * x_D + c3
            D = (x_D, y_D)

            # Validate that the given points lie on the segments
            if (is_on_segment(p1, A, B) and is_on_segment(p2, B, C) and
                is_on_segment(p3, C, D) and is_on_segment(p4, D, A)):
                return [A, B, C, D]
            # Some permutations of vertices might be needed depending on point ordering
            if (is_on_segment(p1, D, A) and is_on_segment(p2, A, B) and
                is_on_segment(p3, B, C) and is_on_segment(p4, C, D)):
                 return [A, B, C, D]

    return None # No solution found

# The four given points
points = [
    (0.3511, 0.2027),
    (0.6753, 0.8303),
    (-0.2845, 0.9905),
    (-0.128, 0.2218)
]

vertices = solve_square_vertices(points)

if vertices:
    # Sort vertices by x-coordinate
    vertices.sort(key=lambda v: v[0])
    
    print("The coordinates of the vertices of the square are:")
    for v in vertices:
        # Print each vertex coordinate pair, formatted to 2 decimal places.
        print(f"({v[0]:.2f}, {v[1]:.2f})")
else:
    print("Could not find a valid square for the given points.")
