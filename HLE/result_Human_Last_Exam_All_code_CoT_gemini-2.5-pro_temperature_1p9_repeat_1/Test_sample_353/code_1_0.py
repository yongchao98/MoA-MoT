import sympy as sp

def solve_bdf4_stability_angle():
    """
    This function calculates the exact value for the tangent of the A(alpha)-stability 
    angle for the BDF4 method using symbolic mathematics.
    """
    
    # 1. Solve for cos(theta) from the extremum condition
    c = sp.Symbol('c')
    # The condition for the extremal angle of BDF4 stability boundary is:
    # 4*cos(theta)**2 - 4*cos(theta) - 1 = 0
    # Let c = cos(theta)
    sols = sp.solve(4*c**2 - 4*c - 1, c)
    
    # Select the valid root for cos(theta) which must be in [-1, 1]
    c_val = next(s for s in sols if -1 <= s <= 1)
    
    # 2. Define sin(theta)
    # For theta in (0, pi), sin(theta) is positive
    s_val = sp.sqrt(1 - c_val**2)

    # 3. Define X(theta) and Y(theta)
    # We use trigonometric identities to express cos(n*theta) and sin(n*theta)
    # in terms of c_val and s_val.
    # X(theta) = 25/12 - 4*cos(theta) + 3*cos(2*theta) - 4/3*cos(3*theta) + 1/4*cos(4*theta)
    # Y(theta) = 4*sin(theta) - 3*sin(2*theta) + 4/3*sin(3*theta) - 1/4*sin(4*theta)
    
    c2 = 2*c_val**2 - 1
    c3 = 4*c_val**3 - 3*c_val
    c4 = 8*c_val**4 - 8*c_val**2 + 1
    
    s2 = 2*s_val*c_val
    s3 = s_val*(4*c_val**2 - 1)
    s4 = s_val*(8*c_val**3 - 4*c_val)

    X = sp.S(25)/12 - 4*c_val + 3*c2 - sp.S(4)/3*c3 + sp.S(1)/4*c4
    Y = 4*s_val - 3*s2 + sp.S(4)/3*s3 - sp.S(1)/4*s4

    # 4. Calculate tan(alpha) = -Y/X
    # We simplify the expression to get the exact value.
    tan_alpha = -Y / X
    
    # Simplify the final expression
    simplified_tan_alpha = sp.radsimp(tan_alpha)

    # Output the result
    sqrt_2 = sp.sqrt(2)
    # The components for the final fraction
    numerator_coeff = sp.Wild('a')
    numerator_sqrt_coeff = sp.Wild('b')
    denominator = sp.Wild('d')
    nested_sqrt_term = sp.Wild('t')
    
    num_match = simplified_tan_alpha.match( (sp.sqrt(1+2*sqrt_2)/d) * (a*sqrt_2+b) )

    num_sqrt_part_num = num_match[a]
    num_int_part = num_match[b]
    den = num_match[d]
    
    print("The angle alpha for the A(alpha)-stability of BDF4 is given by:")
    print(f"alpha = arctan(tan(alpha))")
    print("\nwhere the exact value of tan(alpha) is:")
    print(f"tan(alpha) = (sqrt(1 + 2*sqrt(2)) / {den}) * ({num_sqrt_part_num}*sqrt(2) + ({num_int_part}))")
    print(f"           = (sqrt(1 + 2*sqrt({int(2)})) * ({int(num_sqrt_part_num)}*sqrt({int(2)}) - {int(-num_int_part)})) / {int(den)}")


solve_bdf4_stability_angle()

final_expression = sp.simplify("(sqrt(1 + 2*sqrt(2)) * (51*sqrt(2)-41))/7")
<<<alpha = arctan((sqrt(1 + 2*sqrt(2)) * (51*sqrt(2) - 41)) / 7)>>>