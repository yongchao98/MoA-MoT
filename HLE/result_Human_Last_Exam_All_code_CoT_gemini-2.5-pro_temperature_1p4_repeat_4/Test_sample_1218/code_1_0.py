def solve_k_n_relation():
    """
    Calculates the maximum value of n in terms of k for a k-uniform
    intersecting family with full differences of size k-1.
    """
    # We will use an example value for k, for instance, k=4.
    # The user can change this value to test other cases.
    k = 4
    
    # Based on the derivation, the maximum value of n is given by the formula n = 2k - 1.
    n = 2 * k - 1
    
    print(f"For a given k, the maximum value of n is determined by the equation n = 2k - 1.")
    print(f"For the example case where k = {k}, the calculation is:")
    print(f"n = 2 * {k} - 1 = {n}")

solve_k_n_relation()