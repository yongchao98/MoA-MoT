import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary {prime: exponent}.
    """
    factors = {}
    if n <= 0:
        return factors
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

def solve_for_l(l):
    """
    Calculates and prints the steps for |U_l| and |T_l|.
    """
    print(f"--- For l = {l} ---")
    factors_l = get_prime_factorization(l)
    
    # --- Part A ---
    print("Part A: |U_l|")
    
    l_odd = l
    if 2 in factors_l:
        l_odd = l // (2**factors_l[2])
    
    if l_odd == 1:
        u_l = 2
        print(f"l = {l}. The odd part of l is 1.")
        print(f"|U_l| = 2")
    else:
        l_odd_factors = get_prime_factorization(l_odd)
        u_l = 2
        expr_a = f"|U_l| = 2"
        calc_a = " = 2"
        
        for p, e in l_odd_factors.items():
            term = p**e + p**(e-1)
            u_l *= term
            expr_a += f" * ({p}^{e} + {p}^{e-1})"
            calc_a += f" * ({p**e} + {p**(e-1)})"

        print(f"The odd part of l is {l_odd}.")
        print(f"The formula is |U_l| = 2 * product(p^e + p^(e-1)) for each prime factor p^e of the odd part.")
        print(expr_a)
        print(calc_a)
        print(f" = {u_l}")

    # --- Part B ---
    print("\nPart B: |T_l|")
    
    if l == 1:
        t_l = 1
        print("For l=1, the condition is lambda^2 < max(1, 2), so lambda=1.")
        print("This implies m=1, n=1. T_1 = {(1, 1, 1)}.")
        print(f"|T_l| = 1")
    else:
        print("For l>1, the condition is lambda < l, which excludes the case lambda=l.")
        print("The total number of triples is tau(l^2), so |T_l| = tau(l^2) - 1.")
        
        expr_b = "|T_l| = product(2*e + 1) - 1"
        calc_b = " = ("
        
        t_l = 1
        terms = []
        for p, e in factors_l.items():
            term = 2 * e + 1
            t_l *= term
            terms.append(f"(2*{e} + 1)")
        
        calc_b += " * ".join(terms) + ") - 1"
        
        print(f"The prime factorization of l is {factors_l}.")
        print(expr_b)
        print(calc_b)
        print(f" = {t_l} - 1 = {t_l - 1}")
    print("-" * (13 + len(str(l))))

# Example with a value, e.g., l = 12
l_to_solve = 12
solve_for_l(l_to_solve)

# Example with l=1
solve_for_l(1)

# Example with an odd l > 1
solve_for_l(15)
