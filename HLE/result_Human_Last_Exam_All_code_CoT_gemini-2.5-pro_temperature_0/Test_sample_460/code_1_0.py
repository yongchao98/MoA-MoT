import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve():
    """
    Calculates the smallest integer u based on the given parameters t and m.
    The formula for the smallest u is t * max_{k in {0..m-1}} C(m-1, k).
    """
    t = 20
    m = 4
    
    # We need to find the maximum of C(m-1, k) for k from 0 to m-1.
    n = m - 1
    max_binom_coeff = 0
    for k in range(n + 1):
        coeff = combinations(n, k)
        if coeff > max_binom_coeff:
            max_binom_coeff = coeff
            
    # The smallest u is t times this maximum value.
    u = t * max_binom_coeff
    
    print(f"Given parameters: t = {t}, m = {m}")
    print(f"The formula for the smallest u is t * max(C(m-1, k) for k in [0, m-1]).")
    print(f"For m = {m}, we calculate max(C({m-1}, k)).")
    print(f"The maximum value of C({n}, k) is {max_binom_coeff}.")
    print(f"The final equation is: u = t * {max_binom_coeff}")
    print(f"u = {t} * {max_binom_coeff}")
    print(f"The smallest integer u is: {u}")

solve()