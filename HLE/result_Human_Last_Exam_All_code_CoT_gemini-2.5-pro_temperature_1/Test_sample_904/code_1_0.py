import sympy

def derive_coefficients():
    """
    This function symbolically derives the coefficients A(r) and B(r) for the
    governing linear equation of the fluid interface.
    """
    # 1. Define symbolic variables and functions
    # r: radial position
    # xi(r): displacement of the fluid interface
    # gamma(r): surface tension, which varies with r due to the electric field
    r = sympy.symbols('r', real=True, positive=True)
    xi = sympy.Function('xi')(r)
    gamma = sympy.Function('gamma')(r)

    # 2. State the physical principle
    # The governing equation arises from the linearized force balance on an annular
    # element of the interface. This leads to the generalized Young-Laplace equation:
    # (1/r) * d/dr [ r * gamma(r) * d(xi)/dr ] - Delta_P = 0
    # where Delta_P is the pressure difference (electrostatic + constant).
    # We want to write this in the form: A(r)*xi'' + B(r)*xi' + C(r, xi) = 0

    # 3. Symbolically perform the differentiation to expand the equation.
    # Define the first derivative of xi for clarity
    xi_prime = sympy.diff(xi, r)

    # This is the term inside the main derivative
    term_to_differentiate = r * gamma * xi_prime

    # Differentiate with respect to r using the product rule
    expanded_derivative = sympy.diff(term_to_differentiate, r)

    # The left-hand side of the ODE (without the pressure term) is (1/r) times this
    lhs_operator = (expanded_derivative / r).expand()

    # 4. Identify the coefficients A(r) and B(r) by collecting terms.
    # A(r) is the coefficient of the second derivative of xi, xi''
    # B(r) is the coefficient of the first derivative of xi, xi'
    xi_double_prime = sympy.diff(xi, r, 2)
    
    A_r = lhs_operator.coeff(xi_double_prime)
    B_r = lhs_operator.coeff(xi_prime)

    # 5. Print the final results. The output represents the derived symbolic expressions.
    print("The governing linear equation for the interfacial shape xi(r) is of the form:")
    print("A(r) * xi''(r) + B(r) * xi'(r) + C(r, xi) = 0\n")
    print("Based on the physical derivation, the coefficients are:")
    
    # Use sympy.pretty_print for a more readable mathematical format
    print("\nA(r) =")
    sympy.pretty_print(A_r)
    
    print("\nB(r) =")
    sympy.pretty_print(B_r)

if __name__ == '__main__':
    derive_coefficients()