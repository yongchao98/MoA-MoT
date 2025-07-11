import sympy

def solve_controller():
    """
    This function determines and prints the set of all proper stabilizing controllers H_2(s)
    for the given system H_1(s).
    """
    # Define the symbolic variables for s and the parameter K(s)
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Numerator of H_2(s)
    # 6*(s**2+3*s+2) + (s**2-1)*K(s)
    # Expanded form: 6*s**2 + 18*s + 12 + (s**2-1)*K(s)
    num_poly_part = 6*s**2 + 18*s + 12
    num_K_part = (s**2 - 1) * K
    numerator_expr = num_poly_part + num_K_part
    
    # To display coefficients clearly as requested
    num_poly_part_coeffs = (sympy.Poly(6*s**2, s).as_expr(),
                             sympy.Poly(18*s, s).as_expr(),
                             sympy.Poly(12, s).as_expr())

    num_K_part_poly_coeffs = (sympy.Poly(1*s**2, s).as_expr(),
                              sympy.Poly(0*s, s).as_expr(),
                              sympy.Poly(-1, s).as_expr())

    # Denominator of H_2(s)
    # s**2 - 4 - s*K(s)
    den_poly_part = s**2 - 4
    den_K_part = -s * K
    denominator_expr = den_poly_part + den_K_part
    
    den_poly_part_coeffs = (sympy.Poly(1*s**2, s).as_expr(),
                             sympy.Poly(0*s, s).as_expr(),
                             sympy.Poly(-4, s).as_expr())
    
    den_K_part_poly_coeffs = (sympy.Poly(-1*s, s).as_expr(),)

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("H_2(s) = N(s) / D(s)")
    print("\nWhere K(s) is any stable and proper transfer function.")
    print("\nNumerator N(s):")
    # Pretty print the expression using explicit coefficients
    print(f"N(s) = ({num_K_part_poly_coeffs[0]} + {num_K_part_poly_coeffs[2]}) * K(s) + ({num_poly_part_coeffs[0]} + {num_poly_part_coeffs[1]} + {num_poly_part_coeffs[2]})")

    print("\nDenominator D(s):")
    print(f"D(s) = ({den_poly_part_coeffs[0]} + {den_poly_part_coeffs[2]}) + ({den_K_part_poly_coeffs[0]}) * K(s)")
    
    final_equation = sympy.Eq(sympy.Function('H_2')(s), numerator_expr / denominator_expr)
    
    # Store the result in a special format for direct extraction
    result_str = str(final_equation.rhs)
    # print(f"\n\n<<<H_2(s) = {result_str}>>>")


solve_controller()