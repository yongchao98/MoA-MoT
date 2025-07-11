def solve_markov_chain_problem():
    """
    Searches for the smallest integer r > 1 satisfying the derived equation.
    The equation is g(p) + g(q) - 11 = g(r), where g(x) = 2x^3 + 3x^2 + 6x.
    We are looking for integer solutions (p, q, r) with p, q, r > 1.
    """
    
    def g(x):
        """Calculates the polynomial g(x) = 2x^3 + 3x^2 + 6x."""
        return 2 * x**3 + 3 * x**2 + 6 * x

    # Set a practical upper limit for the search.
    search_limit = 100

    # Pre-compute g(x) values for efficient lookup.
    # The map's limit must be large enough to contain potential q values.
    # A rough analysis shows q can be larger than r, so we extend the map limit.
    g_map_limit = int(search_limit * 1.5)
    g_values_map = {g(i): i for i in range(2, g_map_limit)}
    
    # Iterate through possible values of r starting from 2.
    for r in range(2, search_limit + 1):
        target_sum = g(r) + 11
        
        # Iterate through p. Assume p <= q without loss of generality.
        # This implies g(p) <= g(q), so 2*g(p) <= target_sum.
        p_limit = int((target_sum / 4)**(1/3.0)) + 2 # Heuristic upper bound for p
        for p in range(2, p_limit):
            g_p = g(p)

            if 2 * g_p > target_sum:
                break

            required_g_q = target_sum - g_p
            
            if required_g_q in g_values_map:
                q = g_values_map[required_g_q]
                
                # A solution is found. We print the details and the result.
                print("Yes, such integers exist.")
                print(f"The smallest possible value of r found is {r}.")
                print("\nThis is based on the equation: g(p) + g(q) - 11 = g(r), where g(x) = 2x^3 + 3x^2 + 6x.")
                print("A solution triplet (p, q, r) is:")
                print(f"p = {p}")
                print(f"q = {q}")
                print(f"r = {r}")
                
                print("\nVerification:")
                print(f"g({p}) + g({q}) - 11 = {g_p} + {required_g_q} - 11 = {g_p + required_g_q - 11}")
                print(f"g({r}) = {g(r)}")
                
                print(f"\nFinal Answer: <<< {r} >>>")
                return

    print("No solution found within the search limit.")

solve_markov_chain_problem()