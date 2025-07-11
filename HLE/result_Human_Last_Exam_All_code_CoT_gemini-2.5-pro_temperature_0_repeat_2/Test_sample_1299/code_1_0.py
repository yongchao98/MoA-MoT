import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary of {prime: exponent}.
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

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = {1}
    for p, e in get_prime_factorization(n).items():
        divs.update({d * (p**i) for d in divs for i in range(1, e + 1)})
    return sorted(list(divs))

def calculate_cardinalities(l):
    """
    Calculates and prints the cardinalities of U_l and T_l.
    """
    if not isinstance(l, int) or l <= 0:
        print("Please provide a positive integer for l.")
        return

    print(f"--- Calculations for l = {l} ---")
    
    l_factors = get_prime_factorization(l)
    exponents = list(l_factors.values())
    
    # --- Part B: Calculate |T_l| ---
    print("\nB) Calculating |T_l|:")
    if l == 1:
        card_T_l = 1
        print("|T_1| is defined as 1 based on the condition lambda^2 < max(1, 2).")
    else:
        prod = 1
        prod_str_parts = []
        for e in exponents:
            prod *= (1 + 2 * e)
            prod_str_parts.append(f"(1 + 2*{e})")
        
        card_T_l = prod - 1
        
        print(f"l = {l} has prime factorization with exponents e_i = {exponents}")
        print(f"|T_{l}| = ({' * '.join(prod_str_parts)}) - 1 = {prod} - 1 = {card_T_l}")

    # --- Part A: Calculate |U_l| ---
    print("\nA) Calculating |U_l|:")
    card_U_l = 0
    divisors_of_l = get_divisors(l)
    
    sum_terms = []
    sum_str_parts = []
    
    for lam in divisors_of_l:
        if lam**2 < 2 * l:
            k = l // lam
            k_factors = get_prime_factorization(k)
            omega_k = len(k_factors)
            term = 2**omega_k
            card_U_l += term
            sum_terms.append(term)
            sum_str_parts.append(f"2^omega({l}/{lam})")

    print(f"The condition is lambda^2 < 2*l = {2*l}.")
    print(f"Divisors of {l} satisfying this are: {[lam for lam in divisors_of_l if lam**2 < 2*l]}")
    
    # Build the equation string
    eq_str1 = f"|U_{l}| = {' + '.join(sum_str_parts)}"
    
    omega_values_str = []
    for lam in [d for d in divisors_of_l if d**2 < 2*l]:
         k = l // lam
         omega_k = len(get_prime_factorization(k))
         omega_values_str.append(f"2^{omega_k}")
    
    eq_str2 = f" = {' + '.join(omega_values_str)}"
    eq_str3 = f" = {' + '.join(map(str, sum_terms))}" if len(sum_terms) > 1 else ""
    eq_str4 = f" = {card_U_l}"
    
    print(eq_str1 + eq_str2 + eq_str3 + eq_str4)

# Example usage with a value for l
l_value = 12
calculate_cardinalities(l_value)
