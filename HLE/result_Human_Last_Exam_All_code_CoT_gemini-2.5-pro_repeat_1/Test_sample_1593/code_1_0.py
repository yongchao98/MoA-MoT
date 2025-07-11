def solve_query_complexity():
    """
    This function calculates and prints the query complexity for the two regimes.
    The analysis shows that for both regimes, the optimal strategy yields a complexity
    of Theta(N*log(N)). This is then converted to the specified (a,b,c) format.
    """

    # For a complexity of Theta(N*log(N)), we are looking for a, b, c such that:
    # sqrt(N^a * (log N)^b * (log log N)^c) = N * log(N)
    # Squaring both sides:
    # N^a * (log N)^b * (log log N)^c = N^2 * (log N)^2
    # This gives a=2, b=2, c=0.
    a = 2
    b = 2
    c = 0

    # The analysis concluded that this complexity applies to both regimes.
    regime1_complexity = (a, b, c)
    regime2_complexity = (a, b, c)

    # Format the output string as requested.
    # We output each number in the final equation by constructing the tuple representation.
    result_string = f"({regime1_complexity[0]},{regime1_complexity[1]},{regime1_complexity[2]}),({regime2_complexity[0]},{regime2_complexity[1]},{regime2_complexity[2]})"
    
    print(result_string)

solve_query_complexity()