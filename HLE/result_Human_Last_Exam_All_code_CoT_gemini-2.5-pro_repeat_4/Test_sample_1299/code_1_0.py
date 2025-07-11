import math

def get_prime_factorization(n):
    """
    Returns a dictionary of the prime factorization of n.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
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

def solve_cardinalities(ell):
    """
    Calculates and prints the cardinalities of U_l and T_l for a given l.
    """
    if not isinstance(ell, int) or ell < 1:
        print("Please provide a positive integer for l.")
        return

    # Part A: |U_l| is simply l.
    answer_A = ell

    # Part B: |T_l| calculation.
    # The formula is (Product_{i=1 to s} (1 + 2*e_i)) - 1 + floor(1/l)
    # where l = p_1^e_1 * ... * p_s^e_s.
    
    if ell == 1:
        # For l=1, result is 1.
        answer_B_val = 1
        equation_B = f"l = 1. |T_l| = {answer_B_val}"
    else:
        factors = get_prime_factorization(ell)
        exponents = list(factors.values())

        # Build the equation string for the output
        factor_str_parts = []
        for p, e in sorted(factors.items()):
            factor_str_parts.append(f"{p}^{e}")
        factor_str = " * ".join(factor_str_parts)
        
        term_strs = []
        product_val = 1
        for e in exponents:
            term_val = 2 * e + 1
            term_strs.append(f"({2}*{e}+1)")
            product_val *= term_val
        
        calc_str = "".join(term_strs)
        answer_B_val = product_val - 1
        
        equation_B = f"l = {factor_str}. |T_l| = {calc_str} - 1 = {answer_B_val}"
    
    # Print the final answer in the requested format.
    print(f"For l = {ell}:")
    print(f"A){answer_A} B){equation_B}")

# --- User can change the value of `l` below ---
# Example values:
# l = 1
# l = 12
# l = 7
l = 60

solve_cardinalities(l)
