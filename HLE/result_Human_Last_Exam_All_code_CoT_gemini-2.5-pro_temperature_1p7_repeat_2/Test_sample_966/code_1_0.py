import sympy

def solve_cohomology():
    """
    Calculates the dimension of the middle cohomology group H^100(X, Q) for a
    complete intersection X of two quadrics in CP^102.
    """
    # Parameters of the variety
    # Ambient space is CP^n
    n = 102
    # Degrees of the defining polynomials
    degrees = [2, 2]
    # Dimension of X is m = n - k
    m = n - len(degrees)
    
    # --- Step 1: Calculate the Euler characteristic chi(X) ---
    
    # The degree of X is the product of the degrees of the defining polynomials
    deg_X = 1
    for d in degrees:
        deg_X *= d
    
    # We use a symbolic variable 'g' for the hyperplane class
    g = sympy.Symbol('g')
    
    # The generating function for the Chern classes needed is (1+g)^(n+1) / product(1+d_i*g)
    numerator = (1 + g)**(n + 1)
    denominator = (1 + degrees[0] * g) * (1 + degrees[1] * g)
    gen_func = numerator / denominator
    
    # We need the coefficient of g^m. We perform a series expansion up to O(g^(m+1)).
    series_expansion = sympy.series(gen_func, g, 0, m + 1)
    
    # Extract the coefficient of g^m
    coeff_g_m = series_expansion.coeff(g, m)
    
    # Calculate the Euler characteristic
    chi_X = deg_X * coeff_g_m
    
    # --- Step 2: Calculate the middle Betti number b_m(X) ---
    # For a smooth complete intersection of dimension m, chi(X) = m + b_m(X)
    b_m = chi_X - m
    
    # --- Step 3: Print the final equation with all numbers ---
    print("The dimension of the middle cohomology group H^100(X,Q) is b_100(X).")
    print("The formula is: b_100(X) = chi(X) - 100")
    print("where chi(X) = deg(X) * [g^100] ( (1+g)^103 / (1+2g)^2 )")
    print("")
    print("Calculation:")
    print(f"Degree of X = {degrees[0]} * {degrees[1]} = {deg_X}")
    print(f"The coefficient [g^100] is: {coeff_g_m}")
    print(f"chi(X) = {deg_X} * {coeff_g_m} = {chi_X}")
    print("")
    print("Final equation:")
    print(f"b_100(X) = {chi_X} - {m} = {b_m}")

solve_cohomology()