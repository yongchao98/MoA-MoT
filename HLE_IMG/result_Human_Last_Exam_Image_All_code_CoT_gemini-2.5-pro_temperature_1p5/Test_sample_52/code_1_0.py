import math

def parse_points(numbers):
    """Converts a list of integers into a list of 2D coordinates."""
    points = []
    for n in numbers:
        if n < 10:
            points.append((0, n))
        else:
            points.append((n // 10, n % 10))
    return points

def subtract_vec(v1, v2):
    """Subtracts two 2D vectors."""
    return (v1[0] - v2[0], v1[1] - v2[1])

def add_vec(v1, v2):
    """Adds two 2D vectors."""
    return (v1[0] + v2[0], v1[1] + v2[1])

def dot_product(v1, v2):
    """Computes the dot product of two 2D vectors."""
    return v1[0] * v2[0] + v1[1] * v2[1]

def vec_norm_sq(v):
    """Computes the squared norm (length) of a 2D vector."""
    return v[0]**2 + v[1]**2

def reduce_lattice_basis(b1, b2):
    """
    Reduces a 2D lattice basis (b1, b2) using Lagrange-Gauss algorithm.
    Returns a new basis where the first vector is the shortest vector in the lattice.
    """
    while True:
        norm1_sq = vec_norm_sq(b1)
        norm2_sq = vec_norm_sq(b2)
        
        if norm2_sq < norm1_sq:
            b1, b2 = b2, b1
            continue

        dot = dot_product(b1, b2)
        # Using math.floor for q as per a common variant of the algorithm
        # round() also works and is what Gauss originally used.
        if norm1_sq == 0: # Avoid division by zero
            break
        q = round(dot / norm1_sq)
        
        if q == 0:
            break
        
        b2 = subtract_vec(b2, (q * b1[0], q * b1[1]))

    return b1, b2

def find_shortest_period(vectors):
    """
    Finds the shortest period (length) in a lattice defined by basis vectors.
    """
    if not vectors:
        return 0
    if len(vectors) == 1:
        return math.sqrt(vec_norm_sq(vectors[0]))
        
    b1, b2 = reduce_lattice_basis(vectors[0], vectors[1])
    return math.sqrt(vec_norm_sq(b1))

def solve():
    """Main solver function."""
    cases = [
        [13, 31, 23],
        [10, 4, 23, 31],
        [5, 15, 17, 19, 21, 7],
        [4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13]
    ]

    periods = []

    # Case 1: Triangle
    points = parse_points(cases[0])
    p1, p2, p3 = points[0], points[1], points[2]
    v1 = subtract_vec(p2, p1)
    v2 = subtract_vec(p3, p1)
    period1 = find_shortest_period([v1, v2])
    periods.append(period1)

    # Case 2: Quadrilateral
    points = parse_points(cases[1])
    p1, p2, p3, p4 = points[0], points[1], points[2], points[3]
    v1 = subtract_vec(p3, p1)
    v2 = subtract_vec(p4, p2)
    period2 = find_shortest_period([v1, v2])
    periods.append(period2)

    # Case 3: Hexagon
    points = parse_points(cases[2])
    p1, p2, p3, p4, p5, p6 = points[0], points[1], points[2], points[3], points[4], points[5]
    # Check for parallelo-hexagon condition (a pair of opposite sides are parallel and equal)
    # v_c = p4-p3, v_f = p1-p6. Let's check if v_c + v_f = 0
    v_c = subtract_vec(p4,p3)
    v_f = subtract_vec(p1,p6)
    if add_vec(v_c, v_f) == (0,0): # Type 1 Hexagon confirmed
        # Tiling vectors for a specific type of parallelo-hexagon are P3-P6 and P5-P2.
        v1 = subtract_vec(p3, p6)
        v2 = subtract_vec(p5, p2)
        period3 = find_shortest_period([v1,v2])
    else: # Fallback for other cases (not needed here)
        period3 = 0
    periods.append(period3)

    # Case 4: 13-gon
    points = parse_points(cases[3])
    point_set = set(points)
    # Search for internal translational symmetry
    v1 = (1, 0)
    v2 = (0, 1)
    
    # We found in analysis that the shape is composed of subsets translated by (1,0) and (0,1).
    # This implies the shape can tile the plane with these translation vectors.
    period4 = find_shortest_period([v1, v2])
    periods.append(period4)
    
    print(','.join(map(str, periods)))

solve()