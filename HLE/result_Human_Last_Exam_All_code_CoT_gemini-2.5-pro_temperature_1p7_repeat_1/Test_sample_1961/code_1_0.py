def h_func(k):
    """Calculates the value of the polynomial H(k) = 2k^3 + 3k^2 + 6k."""
    return 2 * k**3 + 3 * k**2 + 6 * k

def find_smallest_r():
    """
    Finds the smallest integer r > 1 that satisfies H(p) + H(q) - 11 = H(r)
    for some integers p, q > 1.
    """
    r_limit = 100 # Search for r up to this limit
    pq_limit = 100 # Search for p, q up to this limit
    
    # Pre-calculate H(r) values for efficient lookup
    h_values = {h_func(r): r for r in range(2, r_limit + 1)}

    for r_candidate in range(2, r_limit + 1):
        target = h_func(r_candidate) + 11
        found = False
        for p in range(2, pq_limit + 1):
            for q in range(p, pq_limit + 1):
                if h_func(p) + h_func(q) == target:
                    # A solution is found
                    print(f"A solution exists for r = {r_candidate}.")
                    print(f"The values are p = {p}, q = {q}, r = {r_candidate}.")
                    print(f"H({p}) + H({q}) - 11 = {h_func(p)} + {h_func(q)} - 11 = {h_func(p) + h_func(q) - 11}")
                    print(f"H({r_candidate}) = {h_func(r_candidate)}")
                    print(f"Final Answer: The smallest possible value of r is {r_candidate}")
                    return r_candidate
    
    print("No solution found within the search limits.")
    return None

smallest_r = find_smallest_r()
if smallest_r is None:
     # This part will not be reached as a solution is found
     print("The search did not yield a solution. The answer may lie beyond the search range, or no solution exists.")
