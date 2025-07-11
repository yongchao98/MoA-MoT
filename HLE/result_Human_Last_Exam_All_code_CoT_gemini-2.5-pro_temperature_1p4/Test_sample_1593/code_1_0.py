def solve():
    """
    Calculates the query complexity for the two specified regimes.
    The complexity is represented by integers (a, b, c) for the class
    Theta(sqrt(N^a * (log N)^b * (log log N)^c)).

    Our analysis shows that for both regimes, the optimal sorting strategy
    results in a query complexity of Theta(N * log N).

    To find the tuple (a, b, c) for this complexity:
    Query Complexity Q = N * log N
    Q^2 = (N * log N)^2 = N^2 * (log N)^2
    We match this to the format: N^a * (log N)^b * (log log N)^c
    This gives:
    a = 2
    b = 2
    c = 0
    """

    # For the regime N = 2^sqrt(L), the complexity is Theta(N*log(N))
    # which corresponds to (a=2, b=2, c=0).
    a1, b1, c1 = 2, 2, 0

    # For the regime N = 2^((log_2 L)^2), the complexity is also Theta(N*log(N))
    # which corresponds to (a=2, b=2, c=0).
    a2, b2, c2 = 2, 2, 0

    # The problem asks for the answer in the format "(a,b,c),(d,e,f)".
    # Each number should be output in the final equation.
    final_answer_string = f"({a1},{b1},{c1}),({a2},{b2},{c2})"
    print(final_answer_string)

solve()
# The final answer in the requested format will be printed by the function call.
# The format required at the end is <<<answer content>>>.
# Based on the code's output, the answer is (2,2,0),(2,2,0).