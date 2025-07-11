import math

def K(x):
    """
    Computes the function K(x) = 2x^3 + 3x^2 + 6x.
    """
    if not isinstance(x, int) or x <= 1:
        return None
    return 2 * x**3 + 3 * x**2 + 6 * x

def solve():
    """
    Searches for integers p, q, r > 1 that satisfy the equation
    K(p) + K(q) - 11 = K(r).
    It searches for the smallest possible value of r.
    """
    # Set a practical search limit for r.
    search_limit = 200
    
    # Pre-compute K(x) values to speed up the search.
    # The dictionary maps K(x) -> x for x in the search range.
    k_values_map = {K(i): i for i in range(2, search_limit)}

    # Iterate on r from 3 upwards to find the smallest solution.
    for r in range(3, search_limit):
        k_r = K(r)
        # For a solution to exist, p and q must be less than r.
        # We iterate q from 2 up to r-1.
        for q in range(2, r):
            k_q = K(q)
            
            # From the equation, we need K(p) = K(r) - K(q) + 11.
            target_k_p = k_r - k_q + 11
            
            # Check if this target K(p) value exists in our pre-computed map.
            if target_k_p in k_values_map:
                p = k_values_map[target_k_p]
                
                # To report each unique solution once, we enforce p <= q.
                if p <= q:
                    # We have found a valid triplet (p, q, r).
                    k_p = K(p)
                    print(f"Found a solution with p={p}, q={q}, r={r}.")
                    print("The derived algebraic equation is: K(p) + K(q) - 11 = K(r)")
                    print("The numbers in this equation are:")
                    print(f"{k_p} + {k_q} - 11 = {k_r}")
                    print(f"Checking: {k_p + k_q - 11} = {k_r}")
                    
                    # Since we are iterating r in increasing order, the first solution
                    # found will give the smallest r.
                    return r
                    
    # If the loops complete without finding a solution, none exists within the limit.
    return "no"

# Run the solver and print the final answer.
result = solve()

if isinstance(result, int):
    print(f"\nThe smallest possible value of r is {result}.")
    print(f"<<<{result}>>>")
else:
    print("\nNo integer solution found for p, q, r within the search limit.")
    print(f"<<<no>>>")
