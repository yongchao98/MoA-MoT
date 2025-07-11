import math

def get_prime_factorization(n):
    """
    Computes the prime factorization of n.
    Returns a dictionary of {prime: exponent}.
    """
    factors = {}
    # Check for factors of 2
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2
    # Check for odd factors
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def solve_for_l(ell):
    """
    Calculates and prints the values for |U_l| and |T_l| for a given l.
    """
    print(f"For l = {ell}:")
    
    # --- Part A ---
    u_ell = ell
    # Keep dividing by 2 until it's odd
    while u_ell > 0 and u_ell % 2 == 0:
        u_ell //= 2
    print(f"A) |U_l| = l_odd = {u_ell}")
    
    # --- Part B ---
    if ell == 1:
        t_ell = 1
        print(f"B) |T_l| = 1 (special case for l=1)")
    else:
        # Get prime factorization l = p_1^e_1 * ... * p_s^e_s
        factors = get_prime_factorization(ell)
        exponents = list(factors.values())
        
        # Calculate tau(l^2) = (2*e_1+1) * ... * (2*e_s+1)
        tau_l_sq = 1
        
        # Build the expression string
        expr_parts = []
        term_values = []
        for e in exponents:
            expr_parts.append(f"(2*{e}+1)")
            term_values.append(2*e+1)
            tau_l_sq *= (2 * e + 1)
        
        t_ell = tau_l_sq - 1
        
        # Build the full expression for T_l
        full_expr_str = "*".join(expr_parts) + " - 1"
        # Build the calculation steps
        calc_step_str = "*".join(map(str, term_values)) + " - 1"
        
        print(f"B) The prime factorization of {ell} is {factors}.")
        print(f"   The exponents are {exponents}.")
        print(f"   The formula for l > 1 is |T_l| = product(2*e_i + 1) - 1.")
        print(f"   |T_l| = {full_expr_str} = {calc_step_str} = {tau_l_sq - 1} - 1 = {t_ell}")

# Demonstrate with a sample value, e.g., l=12
solve_for_l(12)
print("-" * 20)
# Demonstrate with l=1
solve_for_l(1)