def solve_complexity():
    """
    Calculates and prints the complexity tuple (a,b,c) for the two regimes.

    Based on the analysis, the query complexity for both specified regimes simplifies to
    Theta(N * log(N)).

    We need to express this complexity in the form:
    Theta(sqrt(N^a * (log N)^b * (log log N)^c))

    To do this, we square our complexity expression:
    (N * log(N))^2 = N^2 * (log N)^2

    Now we match the exponents with the target form:
    N^a * (log N)^b * (log log N)^c = N^2 * (log N)^2 * (log log N)^0

    This gives us the values for a, b, and c.
    """
    a = 2
    b = 2
    c = 0

    # The format is (a,b,c),(a,b,c), representing the two regimes.
    # Regime 1: N = 2^sqrt(L)
    # Regime 2: N = 2^((log_2 L)^2)
    # Since the complexity is Theta(N*log(N)) for both, the tuple is the same.
    
    regime1_tuple_str = f"({a},{b},{c})"
    regime2_tuple_str = f"({a},{b},{c})"
    
    print(f"{regime1_tuple_str},{regime2_tuple_str}")

solve_complexity()