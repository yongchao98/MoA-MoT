import math

def get_valuations(coefficients):
    """Computes the 2-adic valuation of a list of integer coefficients."""
    valuations = []
    for c in coefficients:
        if c == 0:
            # In theory, valuation is infinity. We skip this for the polygon.
            continue
        v = 0
        if c != 0:
            temp_c = abs(c)
            while temp_c % 2 == 0 and temp_c != 0:
                v += 1
                temp_c //= 2
        valuations.append(v)
    return valuations

def cross_product(p1, p2, p3):
    """
    Calculates the cross product (or turn) for three points.
    Returns:
    > 0 for a counter-clockwise (left) turn
    < 0 for a clockwise (right) turn
    = 0 for collinear points
    """
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])

def get_newton_polygon_vertices(points):
    """
    Computes the lower convex hull (Newton polygon) vertices for a set of points.
    Assumes points are sorted by x-coordinate.
    """
    if len(points) <= 2:
        return points
    
    lower_hull = []
    for p in points:
        while len(lower_hull) >= 2 and cross_product(lower_hull[-2], lower_hull[-1], p) >= 0:
            lower_hull.pop()
        lower_hull.append(p)
    return lower_hull

def solve():
    """
    Solves the problem of finding the number of double points.
    """
    # The polynomial is y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5
    # The coefficients are for powers x^1, x^2, x^3, x^4, x^5
    coeffs = [8, 1, 4, 4, 8]
    
    # Get the 2-adic valuations
    vals = get_valuations(coeffs)
    
    # Create the points (i, v_2(c_i))
    # Exponents are 1-based index + 1
    points = [(i + 1, vals[i]) for i in range(len(vals))]
    
    print("Step 1: The polynomial is f(x) = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x.")
    
    print("\nStep 2: The points (i, v_2(c_i)) for the Newton polygon are:")
    # Sort points by x-coordinate for the convex hull algorithm
    points.sort()
    for p in points:
        print(f"  Exponent {p[0]}: Coefficient valuation is {p[1]}")
    
    # The sorted points are (1,3), (2,0), (3,2), (4,2), (5,3).
    # We pass them to our algorithm to find the lower convex hull vertices.
    
    vertices = get_newton_polygon_vertices(points)
    
    print("\nStep 3: The vertices of the Newton polygon are:")
    for v in vertices:
        print(f"  {v}")
        
    num_vertices = len(vertices)
    # The number of double points is the number of segments, which is num_vertices - 1.
    num_double_points = num_vertices - 1
    
    print(f"\nStep 4: The polygon has {num_vertices} vertices, which form {num_vertices - 1} segments.")
    print("The number of double points in the stable reduction is equal to the number of segments of the Newton polygon.")
    print(f"Final calculation: Number of double points = {num_vertices} - 1 = {num_double_points}")

solve()
<<<2>>>