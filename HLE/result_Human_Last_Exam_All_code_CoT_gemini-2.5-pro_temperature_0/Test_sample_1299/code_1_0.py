import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary {prime: exponent}.
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
       factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def tau_from_factors(factors):
    """
    Calculates the number of divisors (tau function) from prime factors.
    """
    if not factors:
        return 1
    res = 1
    for exp in factors.values():
        res *= (exp + 1)
    return res

def solve_and_explain(l):
    """
    Calculates |U_l| and |T_l| for a given l and explains the steps.
    """
    print(f"--- Solving for l = {l} ---")

    # --- Prime Factorization and Setup ---
    factors = get_prime_factorization(l)
    if l == 1:
        e = []
        p = []
        s = 0
        factor_str = "1"
    else:
        p = sorted(list(factors.keys()))
        e = [factors[pi] for pi in p]
        s = len(p)
        factor_str = " * ".join([f"{pi}^{factors[pi]}" for pi in p])

    print(f"Prime factorization: l = {l} = {factor_str}")
    if l > 1:
        print(f"Exponents: e = {e}")
    
    d = 1 if l % 2 != 0 else 0
    print(f"Parity: l is {'odd' if d==1 else 'even'}, so d = {d}")
    print("-" * 20)

    # --- Part A: Calculation of |U_l| ---
    print("Part A: Calculating |U_l|")
    print("Formula: |U_l| = tau(l) + (1-d) * tau(l/2)")
    
    tau_l = tau_from_factors(factors)
    
    if d == 0:  # l is even
        tau_l_div_2 = tau_from_factors(get_prime_factorization(l // 2))
        U_l = tau_l + tau_l_div_2
        print(f"|U_{l}| = tau({l}) + tau({l//2})")
        print(f"      = {tau_l} + {tau_l_div_2}")
        print(f"      = {U_l}")
    else:  # l is odd
        U_l = tau_l
        print(f"|U_{l}| = tau({l})")
        print(f"      = {tau_l}")
    print("-" * 20)

    # --- Part B: Calculation of |T_l| ---
    print("Part B: Calculating |T_l|")
    if l == 1:
        T_l = 1
        print("For l=1, |T_1| = 1 by direct calculation.")
    else:
        print("Formula (for l>1): |T_l| = (product_{i=1 to s} (1 + 2*e_i)) - 1")
        prod_term = 1
        prod_calc_str_list = []
        prod_val_list = []
        for i in range(s):
            term = 1 + 2 * e[i]
            prod_term *= term
            prod_calc_str_list.append(f"(1 + 2*{e[i]})")
            prod_val_list.append(str(term))
        
        prod_calc_str = " * ".join(prod_calc_str_list)
        prod_val_str = " * ".join(prod_val_list)
        T_l = prod_term - 1

        print(f"|T_{l}| = {prod_calc_str} - 1")
        if s > 1:
            print(f"      = {prod_val_str} - 1")
        print(f"      = {prod_term} - 1")
        print(f"      = {T_l}")
    print("-" * 20)
    print("\n")

if __name__ == '__main__':
    # Example calculation for l = 12
    solve_and_explain(12)
    
    # Example calculation for l = 9 (odd, composite)
    solve_and_explain(9)

    # Example calculation for l = 1 (base case)
    solve_and_explain(1)
