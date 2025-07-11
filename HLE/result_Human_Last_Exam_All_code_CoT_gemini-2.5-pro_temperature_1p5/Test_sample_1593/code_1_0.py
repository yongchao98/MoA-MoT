def solve():
    """
    This function calculates the complexity classes for the two regimes.
    As derived in the explanation, both regimes result in a complexity of Theta(N log N)
    under an optimal radix-sort-based algorithm.

    The complexity Theta(N log N) is represented as Theta(sqrt(N^2 * (log N)^2 * (log log N)^0)),
    which corresponds to the tuple (a=2, b=2, c=0).
    """

    # Complexity for the first regime: N = 2^sqrt(L)
    # The complexity is Theta(N log N), which corresponds to (2,2,0).
    complexity_1_a = 2
    complexity_1_b = 2
    complexity_1_c = 0
    
    # Complexity for the second regime: N = 2^((log_2 L)^2)
    # The complexity is also Theta(N log N) asymptotically, which corresponds to (2,2,0).
    complexity_2_a = 2
    complexity_2_b = 2
    complexity_2_c = 0

    # Format the output string as requested
    output_string = f"({complexity_1_a},{complexity_1_b},{complexity_1_c}),({complexity_2_a},{complexity_2_b},{complexity_2_c})"
    
    print(output_string)

solve()