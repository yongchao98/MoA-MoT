import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents.
    e.g., for n=12, returns {2: 2, 3: 1}
    """
    factors = {}
    if n == 1:
        return {}
    
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

def solve_dessins_problem(l):
    """
    Calculates and prints the solutions for parts A and B for a given integer l.
    """
    if not isinstance(l, int) or l <= 0:
        print("Error: Please provide a positive integer for l.")
        return

    factors = get_prime_factorization(l)
    exponents = list(factors.values())
    s = len(exponents)

    print(f"Solving for l = {l}")
    if l > 1:
        factor_str = " * ".join([f"{p}^{e}" for p, e in factors.items()])
        print(f"The prime factorization is {l} = {factor_str}.")
    
    # --- Part A ---
    # |U_l| = 2^s, where s is the number of distinct prime factors of l.
    val_A = 2**s
    print("\nA) The number of non-isomorphic unicellular regular dessins |U_l|:")
    print(f"   The number of distinct prime factors is s = {s}.")
    print(f"   The formula is |U_l| = 2^s.")
    print(f"   Calculation: |U_{l}| = 2^{s} = {val_A}")

    # --- Part B ---
    # |T_l| = (product of (2*e_i + 1)) - d
    d = 1 if l % 2 != 0 else 0
    
    prod_term = 1
    for e in exponents:
        prod_term *= (2 * e + 1)

    val_B = prod_term - d
    
    print("\nB) The cardinality |T_l|:")
    print(f"   l is {'odd' if d == 1 else 'even'}, so d = {d}.")
    print(f"   The exponents of the prime factors are e_i = {exponents}.")
    print(f"   The formula is |T_l| = (product over i of (2*e_i + 1)) - d.")
    
    if l > 1:
        prod_expr_str = " * ".join([f"(2*{e}+1)" for e in exponents])
        prod_calc_str = " * ".join([str(2*e+1) for e in exponents])
        print(f"   Calculation: |T_{l}| = ({prod_expr_str}) - {d}")
        print(f"                      = ({prod_calc_str}) - {d}")
        print(f"                      = {prod_term} - {d} = {val_B}")
    else: # l = 1
        print(f"   Calculation: |T_{l}| = (empty product) - d")
        print(f"                      = 1 - {d} = {val_B}")

# You can change this value to solve for any positive integer l.
l_value = 60
solve_dessins_problem(l_value)