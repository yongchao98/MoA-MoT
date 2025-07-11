def f(x):
    """Calculates the value of the polynomial 2x^3 + 3x^2 + 6x."""
    return 2 * x**3 + 3 * x**2 + 6 * x

def solve_diophantine_equation():
    """
    Searches for the smallest integer r > 1 that satisfies
    f(p) + f(q) - 11 = f(r) for some integers p, q > 1.
    """
    search_limit = 100  # A reasonable limit for the search space
    
    # Pre-compute f(x) values for efficient lookup
    f_values = {f(x): x for x in range(2, search_limit + 1)}

    # Iterate through possible values of r starting from the smallest
    for r in range(2, search_limit + 1):
        f_r = f(r)
        target_sum_f_p_f_q = f_r + 11
        
        # Iterate through possible values of p
        for p in range(2, search_limit + 1):
            f_p = f(p)
            
            # Since f is increasing, if f(p) alone is already too large, we can stop
            if f_p >= target_sum_f_p_f_q:
                break
            
            # Calculate the required f(q)
            f_q = target_sum_f_p_f_q - f_p
            
            # Check if this f(q) corresponds to an integer q > 1
            if f_q in f_values:
                q = f_values[f_q]
                # We found a valid triplet (p, q, r)
                # Since we iterate r upwards, the first r found is the smallest
                print(f"Found a solution: p = {p}, q = {q}, r = {r}.")
                print(f"The equation is: f({p}) + f({q}) - 11 = {f_p} + {f_q} - 11 = {f_r}")
                print(f"The value of f({r}) is {f_r}, which matches.")
                print(f"The smallest possible value of r is {r}.")
                return

    print("No solution found within the search limit.")

solve_diophantine_equation()
