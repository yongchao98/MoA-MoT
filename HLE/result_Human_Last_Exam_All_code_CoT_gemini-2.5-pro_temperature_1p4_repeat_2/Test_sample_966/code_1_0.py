import math

def calculate_betti_number():
    """
    Calculates the dimension of the middle cohomology group H^100(X, Q) for a
    complete intersection X of two quadrics in CP^102.
    """
    # Parameters from the problem
    n = 102  # Dimension of the ambient projective space CP^n
    degrees = [2, 2]  # Degrees of the defining polynomials
    k = len(degrees) # Number of polynomials

    # Step 1: Calculate the dimension of the variety X
    m = n - k
    
    # We need to find the Betti number b_m(X), where m = 100.
    # The relation is b_m(X) = chi(X) - m, based on properties of such varieties.

    # Step 2: Calculate the Euler characteristic chi(X) using the formula
    # chi(X) = (d1*...*dk) * [z^m] ( (1+z)^(n+1) / product(1+di*z) )
    d_prod = math.prod(degrees)

    # To find the coefficient [z^m], we use a recurrence relation derived from
    # the formula's generating function:
    # C_p = binomial(n+1, p) - 4*C_{p-1} - 4*C_{p-2}.
    # We need to compute C_m.

    C = {}
    # Initialize base cases for p<0
    C[-2] = 0
    C[-1] = 0

    # Iteratively compute C_p up to p = m
    for p in range(m + 1):
        try:
            binom_val = math.comb(n + 1, p)
        except ValueError:
            # math.comb(n, k) raises ValueError if k > n.
            # For p > n+1, the binomial coefficient is 0.
            binom_val = 0
        
        C[p] = binom_val - 4 * C.get(p - 1, 0) - 4 * C.get(p - 2, 0)
    
    coeff_zm = C[m]
    chi_X = d_prod * coeff_zm
    b_m = chi_X - m
    
    print("Step 1: The dimension of the variety X is m = n - k.")
    print(f"m = {n} - {k} = {m}")
    
    print("\nStep 2: The dimension of the middle cohomology group, b_m(X), is related to the Euler characteristic chi(X).")
    print(f"b_{m}(X) = chi(X) - {m}")
    
    print("\nStep 3: The Euler characteristic chi(X) is computed.")
    print(f"chi(X) = ({' * '.join(map(str, degrees))}) * Coeff(z^{m})[ (1+z)^{n+1} / Product(1+d_i z) ]")
    print(f"Coeff(z^{{{m}}}) was calculated using a recurrence relation.")
    
    print("\n--- Calculation Results ---")
    print(f"The coefficient Coeff(z^{{{m}}}) is: {coeff_zm}")
    print(f"The Euler characteristic chi(X) is: {d_prod} * {coeff_zm} = {chi_X}")
    print(f"The dimension of the middle cohomology group is b_{m}(X) = chi(X) - m:")
    print(f"b_{{{m}}}(X) = {chi_X} - {m} = {b_m}")

calculate_betti_number()