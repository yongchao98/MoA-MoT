import sympy

def find_stabilizing_controllers():
    """
    This function calculates the set of all proper stabilizing controllers H2(s)
    for the plant H1(s) = s / (s^2 - 1) using the Youla-Kucera parametrization.
    """
    # Define the symbolic variable 's' and the parameter function K(s)
    s = sympy.symbols('s')
    K = sympy.Function('K')(s)

    # The plant transfer function is H1(s) = s / (s^2 - 1).
    # It has an unstable pole at s = 1.

    # Step 1: Find a coprime factorization H1(s) = N(s) / D(s)
    # We select a stable denominator, (s+1)^2, to create N(s) and D(s)
    # that are stable and proper rational functions.
    N = s / (s + 1)**2
    D = (s**2 - 1) / (s + 1)**2
    
    # We can simplify D(s)
    D = sympy.simplify(D) # D(s) becomes (s-1)/(s+1)

    # Step 2: Solve the Bezout identity N(s)X(s) + D(s)Y(s) = 1
    # This is equivalent to s*X(s) + (s-1)*(s+1)*Y(s) = (s+1)^2
    # A particular solution (X, Y) with stable and proper functions is found to be:
    X = 4
    Y = (s - 1) / (s + 1)
    
    # We can verify this solution:
    # N*X + D*Y = (s/(s+1)**2)*4 + ((s-1)/(s+1))*((s-1)/(s+1))
    #           = (4*s + (s-1)**2) / (s+1)**2
    #           = (4*s + s**2 - 2*s + 1) / (s+1)**2
    #           = (s**2 + 2*s + 1) / (s+1)**2 = 1. Correct.

    # Step 3: Formulate the family of all stabilizing controllers H2(s)
    # H2(s) = (X + D*K) / (Y - N*K)
    
    # Construct the numerator and denominator of H2(s)
    h2_numerator_expr = X + D * K
    h2_denominator_expr = Y - N * K

    # To express H2(s) as a ratio of polynomials in s, we simplify the expression
    # by multiplying the numerator and denominator by a common factor, (s+1)**2.
    final_num = sympy.simplify(h2_numerator_expr * (s+1)**2)
    final_den = sympy.simplify(h2_denominator_expr * (s+1)**2)
    
    # Expand and collect terms with respect to K(s) for clarity
    final_num = sympy.collect(sympy.expand(final_num), K)
    final_den = sympy.collect(sympy.expand(final_den), K)

    # Extract the polynomial parts for formatted printing
    num_poly = sympy.Poly(final_num, K)
    den_poly = sympy.Poly(final_den, K)

    num_part_k = sympy.Poly(num_poly.coeff_monomial(K), s)
    num_part_1 = sympy.Poly(num_poly.coeff_monomial(1), s)
    
    den_part_k = sympy.Poly(den_poly.coeff_monomial(K), s)
    den_part_1 = sympy.Poly(den_poly.coeff_monomial(1), s)

    print("The set of all proper stabilizing controllers H2(s) is given by:")
    print("H2(s) = Num(s) / Den(s)\n")
    print("Where K(s) is any stable and proper transfer function.\n")
    
    # Print the numerator with each coefficient explicitly stated
    print("Numerator Num(s):")
    print(f"({num_part_1.coeff_monomial(s**2)})*s^2 + ({num_part_1.coeff_monomial(s)})*s + ({num_part_1.coeff_monomial(1)}) + (({num_part_k.coeff_monomial(s**2)})*s^2 + ({num_part_k.coeff_monomial(s)})*s + ({num_part_k.coeff_monomial(1)})) * K(s)\n")
    
    # Print the denominator with each coefficient explicitly stated
    print("Denominator Den(s):")
    # Need to handle cases where degree might be lower
    den_c2 = den_part_1.coeff_monomial(s**2) if den_part_1.degree() >= 2 else 0
    den_c1 = den_part_1.coeff_monomial(s) if den_part_1.degree() >= 1 else 0
    den_c0 = den_part_1.coeff_monomial(1)

    den_k_c2 = den_part_k.coeff_monomial(s**2) if den_part_k.degree() >= 2 else 0
    den_k_c1 = den_part_k.coeff_monomial(s) if den_part_k.degree() >= 1 else 0
    den_k_c0 = den_part_k.coeff_monomial(1)
    
    print(f"({den_c2})*s^2 + ({den_c1})*s + ({den_c0}) + (({den_k_c2})*s^2 + ({den_k_c1})*s + ({den_k_c0})) * K(s)")

find_stabilizing_controllers()