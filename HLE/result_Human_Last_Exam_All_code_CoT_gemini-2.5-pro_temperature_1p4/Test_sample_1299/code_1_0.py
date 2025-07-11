import math

def get_prime_factorization(n):
    """Computes the prime factorization of n.
    
    Args:
        n: A positive integer.
        
    Returns:
        A dictionary mapping each prime factor to its exponent.
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def omega(n, prime_factors_cache):
    """Computes the number of distinct prime factors of n.
    
    Args:
        n: A positive integer.
        prime_factors_cache: A dictionary to cache factorization results.
        
    Returns:
        The number of distinct prime factors of n.
    """
    if n in prime_factors_cache:
        return len(prime_factors_cache[n])
    factors = get_prime_factorization(n)
    prime_factors_cache[n] = factors
    return len(factors)

def get_divisors(n):
    """Computes the divisors of n.
    
    Args:
        n: A positive integer.
        
    Returns:
        A sorted list of divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_for_ell(ell):
    """
    Calculates and prints the cardinalities |U_l| and |T_l| for a given integer l.
    The problem assumes a relationship |U_l| = |T_l|.
    """
    print(f"For l = {ell}:")

    # Handle the special case l = 1
    if ell == 1:
        # For l=1, T_1 = {(1,1,1)}, so |T_1|=1.
        # Assuming |U_1| = |T_1|.
        result_str = "A)[1] B)[1]"
        print(result_str)
        print(f'<<<{result_str}>>>')
        return

    # --- Part B: Calculate |T_l| using prime factorization exponents ---
    factors = get_prime_factorization(ell)
    exponents = list(factors.values())
    
    prod_terms_str_list = []
    prod_terms_val_list = []
    prod_val = 1
    for e in sorted(exponents, reverse=True): # sort for consistent output
        prod_terms_str_list.append(f"(2*{e}+1)")
        term_val = 2 * e + 1
        prod_terms_val_list.append(str(term_val))
        prod_val *= term_val
    card_T = prod_val - 1
    
    b_expr_str = "*".join(prod_terms_str_list)
    b_val_str = "*".join(prod_terms_val_list)
    B_str = f"B)[{b_expr_str} - 1 = {b_val_str} - 1 = {prod_val} - 1 = {card_T}]"

    # --- Part A: Calculate |U_l| using sum over divisors ---
    divs = get_divisors(ell)
    prime_factors_cache = {}
    
    omega_vals = [omega(d, prime_factors_cache) for d in divs]
    two_pow_omega_vals = [2**o for o in omega_vals]
    sum_val = sum(two_pow_omega_vals)
    card_U = sum_val - 1

    omega_str_parts = [f"2^{o}" for o in omega_vals]
    omega_str = " + ".join(omega_str_parts)
    
    two_pow_omega_str_parts = [str(t) for t in two_pow_omega_vals]
    two_pow_omega_str = " + ".join(two_pow_omega_str_parts)

    A_str = f"A)[({omega_str}) - 1 = ({two_pow_omega_str}) - 1 = {sum_val} - 1 = {card_U}]"

    final_output = f"{A_str} {B_str}"
    print(final_output)
    print(f'<<<{final_output}>>>')


# As no specific value for l was provided, I will demonstrate the solution for l = 12.
# 12 = 2^2 * 3^1
solve_for_ell(12)