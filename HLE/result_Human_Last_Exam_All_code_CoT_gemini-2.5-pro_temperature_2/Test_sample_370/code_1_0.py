import sympy

def calculate_cross_section_formula():
    """
    Calculates and prints the symbolic formula for the total cross section sigma.
    """
    # Define symbols
    g, E, M, pi = sympy.symbols('g E M pi', real=True, positive=True)

    # Kinematic variables in terms of E and M
    # a corresponds to 2*E**2
    # b corresponds to M**2
    a = 2 * E**2
    b = M**2

    # The result of the definite integral of the squared matrix element over cos(theta) from -1 to 1.
    # The integral is of the form integral( (2 - F(x))^2 ) dx = 8 - 4 * integral(F(x)) + integral(F(x)^2)
    
    # 4 * integral of F(x) from -1 to 1
    integral_F = (8 * b / a) * sympy.log((2 * a + b) / b)
    
    # integral of F(x)^2 from -1 to 1
    integral_F_squared_term1 = (4 * b) / (2 * a + b)
    integral_F_squared_term2 = ((2 * b**2) / (a * (a + b))) * sympy.log((2 * a + b) / b)
    integral_F_squared = integral_F_squared_term1 + integral_F_squared_term2

    # Total value of the kinematic integral part
    kinematic_integral = 8 - integral_F + integral_F_squared
    
    # Prefactor for the total cross section formula
    # The s in the formula is (2E)^2 = 4E^2
    prefactor = g**4 / (128 * pi * (4 * E**2))

    # Total cross section sigma
    sigma = prefactor * kinematic_integral
    
    # Simplify the final expression
    sigma_simplified = sympy.simplify(sigma)
    
    print("The total cross section sigma is:")
    sympy.pprint(sigma_simplified)

    # For verification, we can also display the expanded version which corresponds more directly to the terms derived.
    # Term for 8 in integral
    term1 = g**4 / (64 * pi * E**2)
    # Term from -4*integral(F)
    term2 = - (g**4 / (512 * pi * E**2)) * ((4 * M**2 / E**2) * sympy.log(1 + 4*E**2/M**2))
    # Terms from integral(F^2)
    term3_1 = (g**4 / (512 * pi * E**2)) * (4 * M**2 / (4*E**2 + M**2))
    term3_2 = (g**4 / (512 * pi * E**2)) * (M**4 / (E**2 * (2*E**2+M**2))) * sympy.log(1 + 4*E**2/M**2)
    
    final_eq = term1 + sympy.simplify(term2) + sympy.simplify(term3_1) + sympy.simplify(term3_2)

    # Uncomment the line below to print a more expanded (but less compact) form
    # print("\nExpanded form:")
    # sympy.pprint(final_eq)

calculate_cross_section_formula()
<<<g**4*(2*M**2*(6*M**4 + 20*E**2*M**2 + 16*E**4)*log(1 + 4*E**2/M**2) - 3*M**6 - 30*E**2*M**4 - 96*E**4*M**2 - 64*E**6)/(256*pi*E**4*(-M**2 - 4*E**2)*(M**2 + 2*E**2))>>>