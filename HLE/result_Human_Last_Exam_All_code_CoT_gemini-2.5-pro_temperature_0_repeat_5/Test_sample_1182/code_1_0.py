import math

def get_v2(n):
    """Computes the 2-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    count = 0
    while n % 2 == 0:
        n //= 2
        count += 1
    return count

def get_newton_polygon_slopes(points):
    """
    Computes the slopes and their lengths of the lower convex hull (Newton polygon).
    Points are given as a dictionary {i: v_i}.
    """
    # Convert dict to list of tuples and sort by x-coordinate
    sorted_points = sorted(points.items())
    
    # Find the lower convex hull using Graham scan idea
    hull = []
    for p in sorted_points:
        while len(hull) >= 2 and \
              (hull[-1][1] - hull[-2][1]) * (p[0] - hull[-1][0]) >= \
              (p[1] - hull[-1][1]) * (hull[-1][0] - hull[-2][0]):
            hull.pop()
        hull.append(p)
        
    slopes = []
    for i in range(len(hull) - 1):
        p1 = hull[i]
        p2 = hull[i+1]
        slope_num = p2[1] - p1[1]
        slope_den = p2[0] - p1[0]
        length = slope_den
        slopes.append({'slope': slope_num / slope_den, 'length': length})
        
    return slopes

def solve():
    """
    Calculates the number of double points in the stable reduction of the curve.
    """
    # Step 1: Determine the genus of the curve
    # y^2 = 8x^5 + 4x^4 + 4x^3 + x^2 + 8x
    # The degree of f(x) is 5.
    degree_f = 5
    g = (degree_f - 1) // 2
    
    # Step 2: Analyze the roots of f(x)
    # f(x) = x * (8x^4 + 4x^3 + 4x^2 + x + 8)
    # One root is r_0 = 0.
    # The other four roots are the roots of g(x) = 8x^4 + 4x^3 + 4x^2 + x + 8
    g_coeffs = {4: 8, 3: 4, 2: 4, 1: 1, 0: 8}
    
    # Get 2-adic valuations of coefficients of g(x)
    g_valuations = {i: get_v2(c) for i, c in g_coeffs.items()}
    
    # Get slopes of the Newton polygon for g(x)
    slopes_info = get_newton_polygon_slopes(g_valuations)
    
    # The valuations of the roots of g(x) are the negative of the slopes.
    g_root_valuations = []
    for s in slopes_info:
        for _ in range(s['length']):
            g_root_valuations.append(-s['slope'])
            
    # Full set of root valuations for f(x)
    # v2(r_0) is infinite, but we treat it as a distinct root.
    # Let's denote the roots as r_0, r_1, r_2, r_3, r_4
    # v2(r_0) = inf
    # v2(r_1) = 3.0 (from slope -3)
    # v2(r_2) = v2(r_3) = v2(r_4) = -1.0 (from slope 1)
    
    # Step 3: Identify the root clusters
    # Two roots r_i, r_j are in the same cluster if v2(r_i - r_j) > 0.
    # Cluster 1: {r_0, r_1}. v2(r_0 - r_1) = v2(r_1) = 3 > 0.
    # For r_i where i > 1, v2(r_0 - r_i) = v2(r_i) = -1 <= 0.
    # For r_i where i > 1, v2(r_1 - r_i) = min(v2(r_1), v2(r_i)) = min(3, -1) = -1 <= 0.
    # For r_i, r_j with i,j > 1, we need to check v2(r_i - r_j).
    # These roots have valuation -1. Let r_i = s_i / 2. v2(s_i)=0.
    # The polynomial for s is z^4+z^3+2z^2+z+16=0. Mod 2, this is z(z^3+z^2+1)=0.
    # The roots s_i mod 2 are the roots of z^3+z^2+1=0, which are distinct.
    # So v2(s_i - s_j) = 0 for i != j.
    # Thus, v2(r_i - r_j) = v2((s_i-s_j)/2) = v2(s_i-s_j) - 1 = 0 - 1 = -1 <= 0.
    # So the roots r_2, r_3, r_4 are in separate clusters.
    clusters = [
        {'size': 2, 'name': 'C1={r0, r1}'},
        {'size': 1, 'name': 'C2={r2}'},
        {'size': 1, 'name': 'C3={r3}'},
        {'size': 1, 'name': 'C4={r4}'}
    ]
    N = len(clusters)
    
    # Step 4: Calculate the genus of each component
    component_genera = []
    for c in clusters:
        g_s = (c['size'] - 1) // 2
        component_genera.append(g_s)
    sum_g_s = sum(component_genera)
    
    # Step 5 & 6: Apply the genus formula and solve for delta
    # g = sum(g_S) + delta - N + 1
    # delta = g - sum(g_S) + N - 1
    delta = g - sum_g_s + N - 1
    
    print("The number of double points (nodes) is found using the formula for the genus of the stable reduction:")
    print("g = (sum of component genera) + delta - N + 1")
    print(f"The genus of the curve is g = {g}.")
    print(f"The number of components in the stable reduction is N = {N}.")
    print(f"The genera of the components are {component_genera}, so their sum is {sum_g_s}.")
    print("Plugging these values into the formula:")
    print(f"{g} = {sum_g_s} + delta - {N} + 1")
    print(f"{g} = delta - {N - 1 - sum_g_s}")
    print(f"delta = {g} + {N - 1 - sum_g_s}")
    print(f"So, the number of double points is delta = {delta}.")

solve()
<<<5>>>