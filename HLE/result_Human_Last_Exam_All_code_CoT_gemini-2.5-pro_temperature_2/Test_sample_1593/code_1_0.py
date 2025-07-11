import math

def solve():
    """
    This function calculates and prints the complexity classes for the two regimes.

    As derived in the explanation, both regimes result in a query complexity of Theta(N log N).
    We convert this complexity into the specified (a,b,c) notation.
    Theta(N log N) corresponds to Theta(sqrt(N^2 * (log N)^2 * (log log N)^0)).
    This gives a=2, b=2, c=0 for both regimes.
    """
    # For Regime 1: N = 2^sqrt(L)
    # The complexity is Theta(N * min(sqrt(L), log N))
    # sqrt(L) = log N, so complexity is Theta(N log N) -> (2,2,0)
    complexity_1 = "(2,2,0)"

    # For Regime 2: N = 2^((log_2 L)^2)
    # The complexity is Theta(N * min(sqrt(L), log N))
    # sqrt(L) > log N asymptotically, so complexity is Theta(N log N) -> (2,2,0)
    complexity_2 = "(2,2,0)"

    # The final answer is a comma-separated string of the two results.
    final_answer = f"{complexity_1},{complexity_2}"
    print(final_answer)

solve()