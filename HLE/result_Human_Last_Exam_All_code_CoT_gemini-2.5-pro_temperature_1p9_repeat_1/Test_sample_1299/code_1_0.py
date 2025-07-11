import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary of prime -> exponent.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def solve_for_l(ell):
    """
    Calculates and prints the values for |T_l| and |U_l| with a step-by-step explanation.
    """
    print(f"--- Solving for l = {ell} ---")

    if ell <= 0:
        print("Error: l must be a positive integer.")
        return

    if ell == 1:
        # Based on the definition of T_l, for l=1, we have (m,n,lambda)=(1,1,1).
        # This is because 1=mn\lambda forces m=n=\lambda=1.
        # We check the condition: lambda^2 < max{l^2, 2} => 1^2 < max{1^2, 2} => 1 < 2, which holds.
        # Thus T_1 = {(1,1,1)}, and |T_1|=1.
        # For U_1, since m=n, there is no distinct pair (n,m,lambda). Hence, |U_1|=1.
        T_l = 1
        U_l = 1
        print("For l=1, this is a special case determined directly from the definitions.")
        print("Final Answer for l=1:")
        print(f"A) |U_1| = {U_l}")
        print(f"B) |T_1| = {T_l}")
        return

    # For l > 1
    prime_factors = get_prime_factorization(ell)
    exponents = list(prime_factors.values())

    print(f"The prime factorization of l={ell} is: ", end="")
    factor_str_parts = [f"{p}^{e}" for p, e in prime_factors.items()]
    print(" * ".join(factor_str_parts))
    print(f"The exponents are e_i: {exponents}\n")

    # --- B: Calculate |T_l| ---
    # Formula: |T_l| = (product of (2*e_i + 1)) - 1
    print("B) To find |T_l|, we first compute the product of (2*e_i + 1) for each exponent e_i:")
    prod_terms = []
    prod_val = 1
    for e in exponents:
        term = 2 * e + 1
        prod_terms.append(f"(2*{e}+1)")
        prod_val *= term

    calc_str = f"Product = {' * '.join(prod_terms)}"
    if len(exponents) > 1:
      calc_terms = [str(2 * e + 1) for e in exponents]
      calc_str += f" = {' * '.join(calc_terms)}"
    calc_str += f" = {prod_val}"
    print(calc_str)
    
    T_l = prod_val - 1
    print(f"Then, |T_{ell}| = Product - 1 = {prod_val} - 1 = {T_l}\n")

    # --- A: Calculate |U_l| ---
    # Formula: |U_l| = |T_l| / 2
    # This division is exact because the product is of odd numbers, which is odd. Subtracting 1 makes it even.
    print("A) To find |U_l|, we use the relation |U_l| = |T_l| / 2 for l > 1:")
    U_l = T_l // 2
    print(f"|U_{ell}| = |T_{ell}| / 2 = {T_l} / 2 = {U_l}\n")

    print("Final Answer for l={}:".format(ell))
    print(f"A) |U_{ell}| = {U_l}")
    print(f"B) |T_{ell}| = {T_l}")

# Execute the solver for a demonstration. You can change this value.
solve_for_l(12)
