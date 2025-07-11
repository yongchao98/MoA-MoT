import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents.
    e.g., for n=60, returns {2: 2, 3: 1, 5: 1}
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def solve_problem(l):
    """
    Calculates and prints the values for |U_l| and |T_l|.
    """
    if l < 1:
        print("l must be a positive integer.")
        return

    prime_factors = get_prime_factorization(l)
    s = len(prime_factors)
    exponents = list(prime_factors.values())
    
    # Part A
    U_l = 2**s
    
    # Part B
    prod_term = 1
    for e in exponents:
        prod_term *= (1 + 2 * e)
    
    delta_l_1 = 1 if l == 1 else 0
    T_l = prod_term - 1 + delta_l_1

    # Formatting the expressions for printing
    s_str = str(s)
    e_strs = [str(e) for e in exponents]
    
    print(f"For l = {l}:")
    print(f"Prime factorization exponents: e_i = {e_strs}")
    print(f"Number of distinct prime factors s = {s}")
    
    # Output for A
    print("\nA) |U_l| = 2^s")
    print(f"|U_{l}| = 2^{s_str} = {U_l}")
    
    # Output for B
    prod_expr_str = "".join([f"(1+2*{e})" for e in e_strs])
    if l == 1:
        print("\nB) |T_l| = (product of (1+2e_i)) - 1 + delta_l,1")
        print(f"|T_{l}| = 1 - 1 + 1 = {T_l}")
    else:
        prod_calc_str = "".join([f"({1+2*int(e)})" for e in e_strs])
        print("\nB) |T_l| = (product of (1+2e_i)) - 1 + delta_l,1")
        print(f"|T_{l}| = {prod_expr_str} - 1 + 0 = {prod_calc_str} - 1 = {prod_term} - 1 = {T_l}")


# Example usage with l = 60
l_value = 60
solve_problem(l_value)