import math

def combinations(n, k):
    """
    Computes the binomial coefficient C(n, k).
    """
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_alon_tarsi_k_nn(n):
    """
    Finds the Alon-Tarsi number of K_n,n using Felikson's formula.
    
    The Alon-Tarsi number for a bipartite graph equals its choice number.
    The choice number of K_n,n is k+1, where k is the integer satisfying:
    C(2k-2, k-1) <= n < C(2k, k)
    """
    k = 1
    while True:
        # Lower bound for the inequality range
        lower_binom = combinations(2 * k - 2, k - 1)
        # Upper bound for the inequality range
        upper_binom = combinations(2 * k, k)
        
        if lower_binom <= n < upper_binom:
            print(f"The value of n is {n}.")
            print(f"We are looking for an integer k that satisfies the inequality:")
            print(f"  C(2k-2, k-1) <= n < C(2k, k)")
            print(f"For k = {k}, the inequality is:")
            print(f"  C(2*{k}-2, {k}-1) <= {n} < C(2*{k}, {k})")
            print(f"Which evaluates to:")
            print(f"  {lower_binom} <= {n} < {upper_binom}")
            print(f"This inequality is true.")
            
            alon_tarsi_number = k + 1
            print(f"\nThe value of k is {k}.")
            print(f"The Alon-Tarsi number of K_{{{n},{n}}} is k + 1.")
            print(f"\nFinal Answer: {k} + 1 = {alon_tarsi_number}")
            return alon_tarsi_number
            
        # If n is too large for the current k, increment k and try again
        k += 1
        if k > n + 1: # Safety break to prevent infinite loops on invalid input
            print(f"Could not find a solution for n={n}.")
            return None

# The graph in question is K_1000,1000
n_val = 1000
solve_alon_tarsi_k_nn(n_val)
