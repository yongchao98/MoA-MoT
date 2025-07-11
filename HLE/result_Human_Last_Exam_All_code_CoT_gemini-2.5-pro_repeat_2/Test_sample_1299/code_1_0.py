import math

def get_prime_factorization(n):
    """Computes the prime factorization of n."""
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

def tau(n):
    """Computes the number of divisors of n."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    num_divisors = 1
    for p in factors:
        num_divisors *= (factors[p] + 1)
    return num_divisors

def solve_dessins(l):
    """
    Calculates the cardinalities for |U_l| and |T_l| and prints the steps.
    """
    print(f"Solving for l = {l}")
    factors = get_prime_factorization(l)
    
    # --- Part A Calculation ---
    print("\n--- Part A: Calculating |U_l| ---")
    
    if l % 2 == 0:
        d = 0
        nu2 = factors.get(2, 0)
        m = l // (2**nu2)
    else:
        d = 1
        m = l
        
    m_squared = m**2
    tau_m_squared = tau(m_squared)
    U_l = (2 - d) * tau_m_squared
    
    print(f"The formula for |U_l| is (2-d) * tau(m^2), where m is the odd part of l.")
    print(f"l = {l} is {'even' if d == 0 else 'odd'}, so d = {d}.")
    print(f"The odd part is m = {m}.")
    print(f"|U_{l}| = (2 - {d}) * tau({m}^2) = {2-d} * tau({m_squared})")
    print(f"The number of divisors, tau({m_squared}), is {tau_m_squared}.")
    print(f"So, |U_{l}| = {2-d} * {tau_m_squared} = {U_l}")

    # --- Part B Calculation ---
    print("\n--- Part B: Calculating |T_l| ---")
    
    if l == 1:
        T_l = 1
        print("|T_1| is 1.")
    else:
        exponents = [factors[p] for p in sorted(factors.keys())]
        prod_terms = [2 * e + 1 for e in exponents]
        
        prod_val = 1
        for term in prod_terms:
            prod_val *= term
        T_l = prod_val - 1
        
        l_factors_str = " * ".join([f"{p}^{e}" for p, e in sorted(factors.items())])
        print(f"The formula for |T_l| for l>1 is tau(l^2) - 1.")
        print(f"l = {l} = {l_factors_str}")
        
        prod_calc_str = " * ".join([f"(2*{e}+1)" for e in exponents])
        prod_terms_str = " * ".join([str(t) for t in prod_terms])

        print(f"|T_{l}| = ({prod_calc_str}) - 1")
        print(f"      = ({prod_terms_str}) - 1")
        print(f"      = {prod_val} - 1 = {T_l}")

    print("\n--- Final Answer ---")
    print(f"<<<A){U_l} B){T_l}>>>")

# As no specific value for l was given, I will use l=60 as a representative example.
# 60 = 2^2 * 3^1 * 5^1
solve_dessins(60)