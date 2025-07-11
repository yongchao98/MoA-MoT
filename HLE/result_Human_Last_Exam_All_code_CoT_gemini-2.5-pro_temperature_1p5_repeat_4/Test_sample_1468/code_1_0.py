import math

def calculate_lower_bound(N, q):
    """
    Calculates the asymptotic lower bound for m.
    The bound is Omega(N * min(q, N/q)).
    This function computes the value inside the Omega notation.
    """
    if not (isinstance(N, int) and isinstance(q, int) and N > 0 and q > 0):
        print("N and q must be positive integers.")
        return

    # Problem constraints
    if not (q <= N / 2):
        print(f"Constraint q <= N/2 is not met. q={q}, N={N}")
        # This is not a fatal error for the formula, but we note it.
        # return

    # The dimensionality of the span of difference vectors for a single y_i
    # is min(q, floor(N/q) - 1).
    # Let's use the asymptotic version: min(q, N/q)
    
    min_val = min(q, N / q)
    
    lower_bound_val = N * min_val
    
    # We present the result in asymptotic notation components
    
    n_div_q = N / q
    
    print(f"Given N = {N} and q = {q}:")
    print(f"The number of rows is N = {N}.")
    print(f"The sparsity is q = {q}.")
    print(f"The term N/q is approximately {n_div_q:.2f}.")
    print(f"We compare q and N/q: min({q}, {n_div_q:.2f}) = {min(q, n_div_q):.2f}.")
    
    final_complexity_term = ""
    if q < math.sqrt(N):
        final_complexity_term = f"N*q = {N*q}"
    elif q > math.sqrt(N):
        final_complexity_term = f"N^2/q = {N**2 / q:.0f}"
    else:
        final_complexity_term = f"N*sqrt(N) = {N*math.sqrt(N):.0f}"
        
    print(f"The asymptotic lower bound for m is of the order of N * min(q, N/q).")
    print(f"For the given values, this evaluates to {N} * {min(q, n_div_q):.2f} = {lower_bound_val:.2f}.")
    print(f"In asymptotic notation, this is Omega({final_complexity_term}).")

# Example usage with some values
# Regime q < sqrt(N)
calculate_lower_bound(N=10000, q=50)
print("\n" + "="*30 + "\n")
# Regime q > sqrt(N)
calculate_lower_bound(N=10000, q=400)
print("\n" + "="*30 + "\n")
# Regime q = sqrt(N)
calculate_lower_bound(N=10000, q=100)