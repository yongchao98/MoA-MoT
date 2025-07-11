import math

def count_double_points():
    """
    Calculates the number of double points in the stable reduction of the curve
    y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5 above p=2.
    """
    # Step 1: Identify the properties of the curve.
    # The curve is given by y^2 = f(x), where f(x) = 8*x^5 + 4*x^4 + 4*x^3 + x^2 + 8*x.
    # The degree of f(x) is 5.
    # The genus of a hyperelliptic curve y^2=f(x) is g = floor((deg(f)-1)/2).
    degree = 5
    g = (degree - 1) // 2
    print(f"The curve has genus g = {g}.")
    print("-" * 20)

    # Step 2: Find the branch points.
    # The branch points are the roots of f(x) and, since the degree is odd, the point at infinity.
    # f(x) = x * (8*x^4 + 4*x^3 + 4*x^2 + x + 8).
    # One branch point is x_0 = 0.
    # Let P(x) = 8*x^4 + 4*x^3 + 4*x^2 + x + 8. The other four branch points are the roots of P(x).
    print("The branch points are {0, roots of P(x), infinity}, where P(x) = 8*x^4 + 4*x^3 + 4*x^2 + x + 8.")
    print("-" * 20)

    # Step 3: Determine the 2-adic valuations of the branch points using the Newton polygon.
    # The coefficients of P(x) are a4=8, a3=4, a2=4, a1=1, a0=8.
    # The 2-adic valuations v2(c) are: v2(8)=3, v2(4)=2, v2(1)=0.
    # The points for the Newton polygon are (i, v2(a_i)): (4,3), (3,2), (2,2), (1,0), (0,3).
    print("The points for the Newton polygon of P(x) are (4,3), (3,2), (2,2), (1,0), (0,3).")
    # The lower convex hull of these points determines the valuations of the roots.
    # Segment 1: from (0,3) to (1,0). Slope = (0-3)/(1-0) = -3. Horizontal length = 1.
    # This implies one root, alpha_1, with valuation v2(alpha_1) = -(-3) = 3.
    v_alpha_1 = 3
    print(f"One root alpha_1 has 2-adic valuation {v_alpha_1}.")
    # Segment 2: from (1,0) to (3,2). Slope = (2-0)/(3-1) = 1. Horizontal length = 2.
    # This implies two roots, alpha_2 and alpha_3, with valuation v2 = -(1) = -1.
    v_alpha_23 = -1
    print(f"Two roots alpha_2, alpha_3 have 2-adic valuation {v_alpha_23}.")
    # Segment 3: from (3,2) to (4,3). Slope = (3-2)/(4-3) = 1. Horizontal length = 1.
    # This implies one root, alpha_4, with valuation v2 = -(1) = -1.
    v_alpha_4 = -1
    print(f"One root alpha_4 has 2-adic valuation {v_alpha_4}.")
    print("The set of branch points is B = {0, alpha_1, alpha_2, alpha_3, alpha_4, infinity}.")
    print("-" * 20)

    # Step 4: Cluster the branch points.
    # Two branch points 'a' and 'b' are in the same cluster if v2(a-b) > 0.
    # Cluster S0: v2(0 - alpha_1) = v2(alpha_1) = 3 > 0. So {0, alpha_1} form a cluster.
    S0_size = 2
    print(f"Cluster S0 = {{0, alpha_1}} has size {S0_size}.")
    # Cluster S_inf: {alpha_2, alpha_3, alpha_4, infinity}.
    # In the coordinate t=1/x, these points correspond to {1/alpha_2, 1/alpha_3, 1/alpha_4, 0}.
    # v2(1/alpha_i) = -v2(alpha_i) = 1 for i=2,3,4.
    # Since these t-values all have positive valuation, they reduce to 0 mod 2.
    # Thus v2(t_i - t_j) > 0 for any pair, so they all cluster together.
    S_inf_size = 4
    print(f"Cluster S_inf = {{alpha_2, alpha_3, alpha_4, infinity}} has size {S_inf_size}.")
    # The number of components in the stable reduction is the number of clusters.
    m = 2
    print(f"The stable reduction has {m} components.")
    print("-" * 20)
    
    # Step 5: Calculate the genera of the components.
    # The genus of a component from a cluster S is floor((|S|-1)/2).
    g0 = (S0_size - 1) // 2
    print(f"The genus of the component C0 (from S0) is g0 = ({S0_size}-1)//2 = {g0}.")
    g_inf = (S_inf_size - 1) // 2
    print(f"The genus of the component C_inf (from S_inf) is g_inf = ({S_inf_size}-1)//2 = {g_inf}.")
    print("-" * 20)

    # Step 6: Calculate the number of double points (delta).
    # The genus of the original curve is related to the component genera by:
    # g = (sum of g_i) + delta - m + 1
    # We solve for delta: delta = g - (sum of g_i) - 1 + m
    sum_g_i = g0 + g_inf
    delta = g - sum_g_i + m - 1
    
    print("The number of double points, delta, is given by the formula:")
    print("delta = g - (g0 + g_inf) + m - 1")
    print(f"delta = {g} - ({g0} + {g_inf}) + {m} - 1")
    print(f"delta = {g} - {sum_g_i} + {m-1}")
    print(f"delta = {delta}")
    
count_double_points()