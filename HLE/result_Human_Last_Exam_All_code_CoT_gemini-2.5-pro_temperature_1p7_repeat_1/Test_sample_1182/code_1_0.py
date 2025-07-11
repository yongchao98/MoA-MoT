import math

def get_valuation(n, p):
    """Computes the p-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    if n % p != 0:
        return 0
    return 1 + get_valuation(n // p, p)

def solve():
    """
    Solves for the number of double points in the stable reduction of the curve.
    """
    p = 2
    # The curve is y^2 = f(x) = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x.
    # One branch point is at x=0. The others are roots of g(x) = f(x)/x.
    # g(x) = 8*x^4 + 4*x^3 + 4*x^2 + x + 8
    # Coefficients of g(x) from x^0 to x^4
    coeffs_g = [8, 1, 4, 4, 8]
    
    print("Step 1: Analyze the polynomial g(x) = 8*x^4 + 4*x^3 + 4*x^2 + x + 8.")
    print(f"Coefficients (from constant term to x^4): {coeffs_g}")

    # Get 2-adic valuations of the coefficients
    valuations = [get_valuation(c, p) for c in coeffs_g]
    print(f"2-adic valuations of coefficients: {valuations}")
    
    # Points for Newton Polygon: (i, v_p(a_i))
    points = list(enumerate(valuations))
    print(f"Points for the Newton Polygon: {points}")

    print("\nStep 2: Determine valuations of roots of g(x) from its Newton Polygon.")
    # The Newton Polygon of g(x) for p=2 has two segments:
    # 1. From (0, 3) to (1, 0), with slope -3 and horizontal length 1.
    #    This implies one root, let's call it e1, with v2(e1) = -(-3) = 3.
    # 2. From (1, 0) to (4, 3), with slope 1 and horizontal length 3.
    #    This implies three roots, e2, e3, e4, with v2(ei) = -(1) = -1.
    root_valuations = {'e1': 3, 'e2': -1, 'e3': -1, 'e4': -1}
    print(f"Valuations of roots of g(x): {root_valuations}")
    
    # The finite branch points of the curve are 0 and the roots of g(x).
    print("\nStep 3: Identify clusters of branch points.")
    print("The branch points are {0, e1, e2, e3, e4, infinity}.")
    # We check the 2-adic valuation of the difference between pairs of branch points.
    # A positive valuation indicates they are in the same cluster.
    
    # Check v2(e1 - 0)
    v2_e1_minus_0 = root_valuations['e1']
    print(f"v2(e1 - 0) = v2(e1) = {v2_e1_minus_0}")

    # Check other pairs involving 0
    # For i=2,3,4, v2(ei - 0) = v2(ei) = -1, which is not positive.
    # Check pairs not involving 0 or infinity
    # v2(e1 - ei) for i=2,3,4 = min(v2(e1), v2(ei)) = min(3, -1) = -1.
    # v2(ei - ej) for i,j in {2,3,4} can be shown to be -1 (from the residual polynomial).
    # Check pairs involving infinity. v2(b - infinity) = -v2(b) for finite b.
    # v2(e1 - infinity) = -v2(e1) = -3.
    # v2(ei - infinity) for i=2,3,4 = -v2(ei) = 1.
    
    print("Two branch points b_i, b_j form a cluster if v2(b_i - b_j) > 0.")
    print(" - For the pair {0, e1}, the valuation of the difference is v2(e1 - 0) = 3 > 0.")
    print(" - For pairs {ei, infinity} where i=2,3,4, v2(ei-inf) = -v2(ei) = 1 > 0. These form a cluster.")
    print("These two clusterings are equivalent under a change of coordinates (x -> 1/x). We'll analyze the first one.")
    
    clusters = []
    if v2_e1_minus_0 > 0:
        clusters.append({'points': ('0', 'e1'), 'v_diff': v2_e1_minus_0})

    print(f"\nFound {len(clusters)} non-trivial cluster(s):")
    for cluster in clusters:
        print(f" - Cluster {cluster['points']} with valuation of difference = {cluster['v_diff']}")

    print("\nStep 4: Calculate the number of double points.")
    total_double_points = 0
    print("The number of double points from a 2-point cluster {b_i, b_j} is v2(b_i - b_j) - 1.")
    
    for cluster in clusters:
        contribution = cluster['v_diff'] - 1
        print(f"The cluster {cluster['points']} contributes {cluster['v_diff']} - 1 = {contribution} double point(s).")
        total_double_points += contribution

    print("\nFinal Result:")
    print(f"The total number of double points in the stable reduction is {total_double_points}.")

solve()