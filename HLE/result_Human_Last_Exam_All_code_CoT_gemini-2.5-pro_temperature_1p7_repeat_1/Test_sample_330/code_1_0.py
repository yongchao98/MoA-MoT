import numpy as np

def solve_square_reconstruction():
    """
    Solves the problem of finding the vertices of a square given one point on each side.
    """
    # The coordinates of the four marked points.
    P1 = np.array([0.3511, 0.2027])
    P2 = np.array([0.6753, 0.8303])
    P3 = np.array([-0.2845, 0.9905])
    P4 = np.array([-0.128, 0.2218])

    # There are three distinct ways to pair the points as being on opposite sides.
    # We test each pairing to find the one that forms a square.
    pairings = [
        ((P1, P3), (P2, P4)),  # Pairing 1: P1/P3 opposite, P2/P4 opposite
        ((P1, P2), (P3, P4)),  # Pairing 2: P1/P2 opposite, P3/P4 opposite
        ((P1, P4), (P2, P3))   # Pairing 3: P1/P4 opposite, P2/P3 opposite
    ]

    best_config = None
    min_error = float('inf')

    for p_config in pairings:
        # Pa and Pc lie on one pair of opposite sides.
        # Pb and Pd lie on the other pair.
        Pa, Pc = p_config[0]
        Pb, Pd = p_config[1]

        # V_ac is the vector from Pa to Pc.
        # V_bd is the vector from Pb to Pd.
        V_ac = Pc - Pa
        V_bd = Pd - Pb

        # For a given pairing, there are two possible orientations (slopes) for the square's sides.
        # We solve for the slope 'm' of one pair of sides (e.g., those containing Pb and Pd).
        # The other pair of sides will have a slope of '-1/m'.
        
        # Candidate slope 1
        denominator1 = V_ac[0] - V_bd[1]
        if abs(denominator1) > 1e-9:
            m = (V_ac[1] + V_bd[0]) / denominator1
            m_perp = -1.0 / m
            
            # Determine the constant 'c' for each of the four lines (y = mx + c)
            c_a = Pa[1] - m_perp * Pa[0]
            c_c = Pc[1] - m_perp * Pc[0]
            c_b = Pb[1] - m * Pb[0]
            c_d = Pd[1] - m * Pd[0]

            # The condition for a square is that the side lengths must be equal.
            # side_length = |c1 - c2| / sqrt(m^2 + 1)
            # This simplifies to |c_a - c_c| * |m| = |c_b - c_d|
            error = abs(abs(c_a - c_c) * abs(m) - abs(c_b - c_d))

            if error < min_error:
                min_error = error
                best_config = {'m1': m, 'm2': m_perp, 'points': (Pa, Pb, Pc, Pd)}

        # Candidate slope 2
        denominator2 = V_ac[0] + V_bd[1]
        if abs(denominator2) > 1e-9:
            m = (V_ac[1] - V_bd[0]) / denominator2
            m_perp = -1.0 / m

            c_a = Pa[1] - m_perp * Pa[0]
            c_c = Pc[1] - m_perp * Pc[0]
            c_b = Pb[1] - m * Pb[0]
            c_d = Pd[1] - m * Pd[0]
            
            error = abs(abs(c_a - c_c) * abs(m) - abs(c_b - c_d))

            if error < min_error:
                min_error = error
                best_config = {'m1': m, 'm2': m_perp, 'points': (Pa, Pb, Pc, Pd)}
    
    # After checking all possibilities, we use the best configuration found.
    m1 = best_config['m1']  # Slope for sides through Pb, Pd
    m2 = best_config['m2']  # Slope for sides through Pa, Pc
    Pa, Pb, Pc, Pd = best_config['points']

    # Define the four lines of the square's sides.
    # The line through Pa has slope m2.
    c_a = Pa[1] - m2 * Pa[0]
    # The line through Pb has slope m1.
    c_b = Pb[1] - m1 * Pb[0]
    # The line through Pc has slope m2.
    c_c = Pc[1] - m2 * Pc[0]
    # The line through Pd has slope m1.
    c_d = Pd[1] - m1 * Pd[0]

    def find_intersection(m_i, c_i, m_j, c_j):
        # Solves for the intersection of two lines: y = mi*x + ci and y = mj*x + cj
        x = (c_j - c_i) / (m_i - m_j)
        y = m_i * x + c_i
        return np.array([x, y])

    # The vertices are the intersections of adjacent side-lines.
    V1 = find_intersection(m2, c_a, m1, c_b)
    V2 = find_intersection(m1, c_b, m2, c_c)
    V3 = find_intersection(m2, c_c, m1, c_d)
    V4 = find_intersection(m1, c_d, m2, c_a)
    
    vertices = [V1, V2, V3, V4]
    
    # Sort vertices based on their x-coordinate.
    vertices.sort(key=lambda v: v[0])
    
    # Print the final coordinates formatted to two decimal places.
    print("The coordinates of the vertices of the square, sorted by x-coordinate, are:")
    for v in vertices:
        print(f"({v[0]:.2f}, {v[1]:.2f})")

solve_square_reconstruction()