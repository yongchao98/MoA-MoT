import math

def p_adic_val(n_in, p):
    """
    Computes the p-adic valuation of an integer.
    """
    n = int(n_in)
    if n == 0:
        return float('inf')
    power = 0
    # We handle the case of n not being an integer, though coeffs are.
    if n != n_in:
        return -1 # Or some other indicator of non-integer. For this problem this won't happen.
    
    # We only care about powers of p in the factorization of n
    if n % p != 0:
        return 0

    while n % p == 0:
        power += 1
        n //= p
    return power

def get_newton_polygon(coeffs, p):
    """
    Calculates the lower convex hull (Newton Polygon) for a polynomial.
    coeffs: a dictionary of degree -> coefficient.
    p: the prime number for the valuation.
    """
    points = []
    for power, coeff in coeffs.items():
        if coeff != 0:
            val = p_adic_val(coeff, p)
            points.append((power, val))

    points.sort()
    
    hull = []
    # Using Andrew's monotone chain algorithm to find the lower convex hull
    for point in points:
        while len(hull) >= 2 and \
              (hull[-1][1] - hull[-2][1]) * (point[0] - hull[-1][0]) >= \
              (point[1] - hull[-1][1]) * (hull[-1][0] - hull[-2][0]):
            hull.pop()
        hull.append(point)
        
    return hull

def analyze_polynomial():
    """
    Analyzes the polynomial from the problem to determine the stable reduction properties.
    """
    # The curve is z^2 = 2*x^5 + 2*x^3 + 1.
    # To find the valuations of the branch points (roots of f(x)=0),
    # we analyze the reciprocal polynomial g(w) = w^5 + 2*w^2 + 2 = 0, where w=1/x.
    # The coefficients are for powers 5, 2, 0.
    g_coeffs = {5: 1, 2: 2, 0: 2}
    prime = 2

    hull = get_newton_polygon(g_coeffs, prime)
    
    print("The curve is z^2 = 2*x^5 + 2*x^3 + 1.")
    print("We analyze the roots of the reciprocal polynomial g(w) = w^5 + 2*w^2 + 2 = 0.")
    print(f"The vertices of the Newton polygon for g(w) with respect to p=2 are: {hull}")

    print("\nThe slopes of the Newton polygon determine the valuations of the roots:")
    for i in range(len(hull) - 1):
        p1 = hull[i]
        p2 = hull[i+1]
        slope_num = p2[1] - p1[1]
        slope_den = p2[0] - p1[0]
        length = p2[0] - p1[0]
        
        # Root valuations are the negative of the slopes.
        root_val_num = -slope_num
        root_val_den = slope_den
        
        print(f"Segment from {p1} to {p2}:")
        print(f"  - Slope: {slope_num}/{slope_den}")
        print(f"  - This implies there are {length} roots w with 2-adic valuation v(w) = {root_val_num}/{root_val_den}.")

    print("\nThis partitions the 5 finite branch points x_i=1/w_i into two clusters based on their valuations.")
    print("A change of variables z = 1+2Z leads to a new model of the curve:")
    print("  (1+2Z)^2 = 2*x^5+2*x^3+1")
    print("  1+4Z+4Z^2 = 2*x^5+2*x^3+1")
    print("  4Z(Z+1) = 2x^3(x^2+1)")
    print("  2*Z*(Z+1) = x^3*(x^2+1)")
    print("\nReducing this equation modulo 2, we get the special fiber:")
    print("  0 = x^3*(x^2+1) mod 2")
    print("  x^3*(x+1)^2 = 0 mod 2")
    print("The special fiber consists of two components: x=0 with multiplicity 3, and x=1 with multiplicity 2.")
    print("The double point of the stable reduction is formed where these components meet after resolution.")
    print("The thickness of this double point is 3.")

analyze_polynomial()