import math

def solve():
    """
    Calculates and prints the complexity tuples for the two regimes.

    My reasoning, as outlined in the plan, leads to the following conclusion:

    For both regimes, N=2^sqrt(L) and N=2^((log L)^2), the optimal query complexity
    is found to be Theta(N log N). This is because simpler strategies like C-sort
    outperform more complex strategies like H-sort or block-based sorts under these
    specific N-L relationships. The costs of these other strategies grow much faster.

    To convert Theta(N log N) to the (a,b,c) notation, we solve:
    sqrt(N^a * (log N)^b * (log log N)^c) = N * log N

    Squaring both sides gives:
    N^a * (log N)^b * (log log N)^c = (N * log N)^2 = N^2 * (log N)^2

    By comparing the exponents of N, log N, and log log N, we find:
    a = 2
    b = 2
    c = 0

    Therefore, the complexity class for both regimes is represented by the tuple (2,2,0).
    """

    # Complexity for the first regime: N = 2^sqrt(L)
    # The analysis shows the optimal complexity is Theta(N log N).
    a1, b1, c1 = 2, 2, 0

    # Complexity for the second regime: N = 2^((log L)^2)
    # The analysis shows the optimal complexity is also Theta(N log N).
    a2, b2, c2 = 2, 2, 0

    # Format the output as specified
    final_answer_string = f"({a1},{b1},{c1}),({a2},{b2},{c2})"
    print(final_answer_string)


solve()