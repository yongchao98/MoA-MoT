def solve_complexity():
    """
    This function analyzes the query complexity for sorting N bitstrings of length L
    in two different regimes and prints the result in the specified (a,b,c) format.

    The analysis shows that for both regimes, the optimal worst-case query
    complexity Q is Theta(N * log(N)).

    We need to find integers (a, b, c) such that:
    Q = Theta(sqrt(N^a * (log N)^b * (log log N)^c))

    Squaring our complexity Q gives:
    Q^2 = Theta((N * log N)^2) = Theta(N^2 * (log N)^2 * (log log N)^0)

    By comparing the exponents with the target form, we get:
    a = 2
    b = 2
    c = 0
    """

    # For Regime 1 (N = 2^sqrt(L)), the complexity is Theta(N*log(N)).
    # This corresponds to (a,b,c) = (2,2,0).
    a1 = 2
    b1 = 2
    c1 = 0

    # For Regime 2 (N = 2^((log L)^2)), the complexity is also Theta(N*log(N)).
    # This corresponds to (a,b,c) = (2,2,0).
    a2 = 2
    b2 = 2
    c2 = 0

    # Format the final answer string as "(a1,b1,c1),(a2,b2,c2)".
    # The problem states to output each number in the final equation.
    # We will print the formatted tuple string.
    result_string = f"({a1},{b1},{c1}),({a2},{b2},{c2})"
    print(result_string)

solve_complexity()