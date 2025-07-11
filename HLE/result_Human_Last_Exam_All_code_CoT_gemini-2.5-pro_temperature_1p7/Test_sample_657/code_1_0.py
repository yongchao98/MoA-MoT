import sympy

def find_stabilizing_controllers():
    """
    This function calculates the set of all proper stabilizing controllers H_2(s)
    for the plant H_1(s) = s / (s^2 - 1) using the Youla-Kucera parametrization.
    """

    # Step 1: Define symbolic variables and the plant H_1(s).
    # K(s) is the free parameter for the parametrization.
    s = sympy.Symbol('s')
    K = sympy.Symbol('K(s)')

    H1 = s / (s**2 - 1)
    
    # The plant has an unstable pole at s=1. The stable pole is at s=-1.

    # Step 2: Factorize the plant H_1(s) into N(s)/D(s) where N(s) and D(s)
    # are stable, proper, coprime rational functions (elements of RH_inf).
    # We use an all-pass filter B(s) to separate the unstable part.
    # Let B(s) = (s-1)/(s+1) corresponding to the unstable pole at s=1.
    # We define N(s) and D(s) such that H_1(s) = N(s)/D(s).
    # N(s) = H_1(s) / B(s)^-1 = H_1(s) * B(s) leads to a wrong result, the definition is H_1 = N * D^-1.
    # A correct factorization is:
    # N(s) = s / (s + 1)**2
    # D(s) = (s - 1) / (s + 1)
    # Check: N/D = [s/(s+1)^2] / [(s-1)/(s+1)] = s/((s+1)(s-1)) = s/(s^2-1) = H_1(s).
    # Both N(s) and D(s) are stable and proper.
    
    N_s = s / (s + 1)**2
    D_s = (s - 1) / (s + 1)

    # Step 3: Solve the Bezout Identity N(s)X(s) + D(s)Y(s) = 1 for stable, proper X, Y.
    # [s/(s+1)^2]*X(s) + [(s-1)/(s+1)]*Y(s) = 1
    # Multiplying by (s+1)^2 gives:
    # s*X(s) + (s^2-1)*Y(s) = (s+1)^2 = s^2 + 2*s + 1
    # By inspection, a simple solution with polynomials (which are a subset of rational functions) is:
    # Y(s) = 1
    # s*X(s) + s^2-1 = s^2 + 2*s + 1  => s*X(s) = 2*s + 2 => X(s) = 2.
    # X(s)=2 and Y(s)=1 are stable and proper.
    
    X_s = 2
    Y_s = 1

    # Step 4: Construct the controller H_2(s) using the Youla-Kucera formula:
    # H_2(s) = (X(s) + D(s)K(s)) / (Y(s) - N(s)K(s))
    
    Numerator_H2 = X_s + D_s * K
    Denominator_H2 = Y_s - N_s * K
    
    # Step 5: Simplify the expression into a single rational function.
    # We use sympy to perform the symbolic algebra.
    # Numerator of H2 becomes: 2 + ((s-1)/(s+1)) * K
    # Denominator of H2 becomes: 1 - (s/((s+1)**2)) * K
    # H2 = [ (2*(s+1) + (s-1)*K)/(s+1) ] / [ ((s+1)**2 - s*K)/((s+1)**2) ]
    # H2 = (s+1) * (2*(s+1) + (s-1)*K) / ((s+1)**2 - s*K)

    num_expr = (s+1) * (2*(s+1) + (s-1)*K)
    den_expr = (s+1)**2 - s*K
    
    # Expand the numerator and denominator to get the polynomial form.
    num_expanded = sympy.expand(num_expr)
    den_expanded = sympy.expand(den_expr)

    # Collect terms with powers of s for a cleaner representation.
    num_final = sympy.collect(num_expanded, s)
    den_final = sympy.collect(den_expanded, s)
    
    # Print the result.
    print("The set of all proper stabilizing controllers H_2(s) for H_1(s) = s/(s^2-1) is given by:")
    print("H_2(s) = Num(s, K(s)) / Den(s, K(s))\n")
    print("where K(s) is any stable and proper rational function.\n")
    
    # Extract coefficients and print the equation with each number.
    num_coeff_s2 = sympy.Poly(num_final, s).coeff_monomial(s**2)
    num_coeff_s1 = sympy.Poly(num_final, s).coeff_monomial(s)
    num_coeff_s0 = sympy.Poly(num_final, s).coeff_monomial(1)
    
    den_coeff_s2 = sympy.Poly(den_final, s).coeff_monomial(s**2)
    den_coeff_s1 = sympy.Poly(den_final, s).coeff_monomial(s)
    den_coeff_s0 = sympy.Poly(den_final, s).coeff_monomial(1)

    print(f"Num(s, K(s)) = ({num_coeff_s2})*s^2 + ({num_coeff_s1})*s + ({num_coeff_s0})")
    print(f"Den(s, K(s)) = ({den_coeff_s2})*s^2 + ({den_coeff_s1})*s + ({den_coeff_s0})")

# Execute the function to get the answer
find_stabilizing_controllers()