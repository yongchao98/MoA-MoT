def solve():
    """
    Calculates the query complexity for the two regimes.
    The analysis shows that in both cases, the complexity is Theta(N log N).
    The complexity class is represented by (a,b,c) where the complexity is
    Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    For Theta(N log N), we have:
    (N log N)^2 = N^2 * (log N)^2
    So, a=2, b=2, c=0.
    """

    # For the first regime, N = 2^sqrt(L), which means L = (log2(N))^2.
    # This corresponds to the complexity Theta(N log N), which is (2,2,0).
    result1 = "(2,2,0)"

    # For the second regime, N = 2^((log2 L)^2), which means L = 2^sqrt(log2(N)).
    # In this regime, L > (log2 N)^2, which means the complexity is also Theta(N log N).
    # This again corresponds to (2,2,0).
    result2 = "(2,2,0)"

    # The final answer should be in the format "(a1,b1,c1),(a2,b2,c2)".
    final_answer = f"{result1},{result2}"
    print(final_answer)

solve()