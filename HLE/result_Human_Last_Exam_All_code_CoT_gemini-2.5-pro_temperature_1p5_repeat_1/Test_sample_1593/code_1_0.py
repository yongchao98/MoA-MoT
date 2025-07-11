def solve():
    """
    This function prints the query complexity for the two regimes.
    Our analysis shows that for both regimes, the optimal algorithm has a query complexity of Theta(N log N).
    
    The query complexity is given in the format (a,b,c) representing the class
    Theta(sqrt(N^a * (log N)^b * (log log N)^c)).

    If the query complexity Q = Theta(N log N), then
    Q^2 = Theta((N log N)^2) = Theta(N^2 * (log N)^2 * (log log N)^0).

    Comparing this with N^a * (log N)^b * (log log N)^c, we get:
    a = 2
    b = 2
    c = 0
    
    This gives the tuple (2,2,0) for both regimes.
    """
    
    # Complexity for the first regime: N = 2^sqrt(L)
    # This corresponds to L = (log N)^2, the crossover point.
    # Complexity is Theta(N log N), which gives the tuple (2,2,0).
    regime1_tuple = "(2,2,0)"
    
    # Complexity for the second regime: N = 2^((log L)^2)
    # For large L, this corresponds to L > (log N)^2.
    # Complexity is Theta(N log N), which gives the tuple (2,2,0).
    regime2_tuple = "(2,2,0)"
    
    final_answer = f"{regime1_tuple},{regime2_tuple}"
    
    print(final_answer)

solve()