def solve():
    """
    Calculates and prints the minimum number of trials T(n) for given values of n.
    The problem asks to output the full equation for each calculation and then
    the final comma-separated list of values.
    """
    values_of_n = [2, 3, 1234, 6712]
    results = []

    for n in values_of_n:
        # The formula for the minimum trials is T(n) = 2n - 1.
        result = 2 * n - 1
        
        # Outputting each number in the final equation as requested.
        print(f"T({n}) = 2 * {n} - 1 = {result}")
        
        results.append(str(result))
    
    # The final deliverable is the set of values, separated by a comma.
    # While the above prints show the calculation, the final answer string is below.
    # To avoid confusion, I will add a clarifying header.
    print("\nThe comma-separated values are:")
    print(",".join(results))

solve()