import math

def p_adic_valuation(n, p):
    """Calculates the p-adic valuation of n, denoted as nu_p(n)."""
    if n == 0:
        return float('inf')
    count = 0
    while n % p == 0 and n != 0:
        count += 1
        n //= p
    return count

def p_adic_valuation_factorial(n, p):
    """Calculates the p-adic valuation of n! using Legendre's formula."""
    if n < 0:
        return float('nan')
    if n == 0:
        return 0
    
    count = 0
    power_of_p = p
    while power_of_p <= n:
        count += n // power_of_p
        power_of_p *= p
    return count

def solve_k_theory_problem():
    """
    Finds the largest natural number n such that the (2n)th K-group of Z/27 is nonzero.
    """
    p = 3
    k = 3
    
    largest_n = 0
    # The exponent E will become negative for large n, so we can limit the search.
    # A search up to n=100 is more than sufficient.
    print("Searching for the largest natural number n...")
    print("n | n+1 | (n/2)! | nu_3(n+1) | nu_3((n/2)!) | E = 2 + nu_3(n+1) - nu_3((n/2)!) | K_{2n}(Z/27) is non-zero?")
    print("-" * 100)

    for n in range(2, 101, 2):
        n_div_2 = n // 2
        
        nu_3_n_plus_1 = p_adic_valuation(n + 1, p)
        nu_3_n_div_2_factorial = p_adic_valuation_factorial(n_div_2, p)
        
        exponent = (k - 1) + nu_3_n_plus_1 - nu_3_n_div_2_factorial
        is_nonzero = exponent > 0
        
        if is_nonzero:
            largest_n = n
        
        print(f"{n:2} | {n+1:3} | {n_div_2:2}!   | {nu_3_n_plus_1:9} | {nu_3_n_div_2_factorial:12} | {exponent:33} | {is_nonzero}")

    print("\n" + "="*40)
    print(f"The largest natural number n is: {largest_n}")
    print("="*40)

solve_k_theory_problem()
