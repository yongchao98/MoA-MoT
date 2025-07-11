import sympy

def solve_fluid_equation():
    """
    This function symbolically derives the coefficients A(r) and B(r)
    for the governing equation of the fluid interface.
    """
    # Define symbols for the physical quantities
    r = sympy.Symbol('r', positive=True) # radial position
    gamma = sympy.Symbol('gamma') # surface tension
    xi = sympy.Function('xi')(r) # interfacial displacement
    P_es = sympy.Function('P_es')(r) # electrostatic pressure

    # From the derivation, the coefficients of the differential equation are identified.
    # The equation is: A(r) * xi''(r) + B(r) * xi'(r) + C(r, xi) = 0

    # Coefficient of the second derivative term
    A_r = gamma

    # Coefficient of the first derivative term
    B_r = gamma / r

    # The term C(r, xi) corresponds to the external pressure, P_es(r).
    C_r_xi = P_es

    # Construct the full governing equation using the derived coefficients
    governing_equation = A_r * xi.diff(r, 2) + B_r * xi.diff(r, 1) + C_r_xi

    # Print the results in a structured way
    print("The governing linear differential equation for the interfacial shape xi(r) is:")
    # Using sympy.pretty_print for a clearer mathematical output
    sympy.pprint(sympy.Eq(governing_equation, 0), use_unicode=True)
    
    print("\n" + "="*50 + "\n")
    print("Comparing this equation to the standard form:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
    
    print("\nWe find the coefficients A(r) and B(r) to be:")

    # Print each coefficient as requested
    print("\nCoefficient A(r):")
    sympy.pprint(sympy.Eq(sympy.Symbol('A(r)'), A_r), use_unicode=True)
    
    print("\nCoefficient B(r):")
    sympy.pprint(sympy.Eq(sympy.Symbol('B(r)'), B_r), use_unicode=True)

if __name__ == "__main__":
    solve_fluid_equation()
