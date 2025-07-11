import math

def get_prime_factorization(n):
    """
    Returns a dictionary where keys are prime factors and values are their exponents.
    e.g., for n=60, returns {2: 2, 3: 1, 5: 1}
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def calculate_T_l_cardinality(l):
    """
    Calculates the cardinality of the set T_l and prints the calculation.
    """
    if not isinstance(l, int) or l < 1:
        print("Error: l must be a positive integer.")
        return

    print(f"Calculating |T_{l}| for l = {l}")
    
    if l == 1:
        print("For l = 1, |T_1| = 1")
        return

    # For l > 1, the formula is |S_l| - 1
    # First, get the prime factorization for l
    prime_factors = get_prime_factorization(l)
    
    # Calculate |S_l| = product of (1 + 2*e_i)
    # And build the expression string
    s_l_val = 1
    expression_parts = []
    
    # Sort primes for consistent output
    sorted_primes = sorted(prime_factors.keys())
    
    for p in sorted_primes:
        e = prime_factors[p]
        term_val = 1 + 2 * e
        s_l_val *= term_val
        expression_parts.append(f"(1 + 2*{e})")
        
    s_l_expression = " * ".join(expression_parts)
    
    # The final cardinality |T_l| is |S_l| - 1
    t_l_val = s_l_val - 1
    
    print(f"The prime factorization of {l} is " + " * ".join([f"{p}^{e}" for p, e in prime_factors.items()]))
    print(f"For l > 1, |T_{l}| = (|S_{l}|) - 1")
    print(f"|S_{l}| = {s_l_expression} = {s_l_val}")
    print(f"So, |T_{l}| = {s_l_val} - 1 = {t_l_val}")

# Example usage with l = 60
calculate_T_l_cardinality(60)