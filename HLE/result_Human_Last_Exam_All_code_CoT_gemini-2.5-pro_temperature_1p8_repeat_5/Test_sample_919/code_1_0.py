import sympy

def display_force_equation():
    """
    This function uses the sympy library to symbolically represent and print
    the final equation for the force per unit area on the conducting plane.
    """
    # Define the symbols used in the equation.
    # We treat them as real and positive, which is physically appropriate.
    mu_0, mu, K_0, a, d, y = sympy.symbols('mu_0 mu K_0 a d y', real=True, positive=True)

    # Build the numerator of the force equation
    # It includes the magnetic constant (mu_0), the current amplitude squared (K_0^2),
    # and the spatial variation term (sin(a*y)^2). The factor of 1/2 is also here.
    numerator_term = (sympy.Rational(1, 2)) * mu_0 * K_0**2 * sympy.sin(a*y)**2

    # Build the denominator of the force equation
    # This part comes from applying boundary conditions and includes hyperbolic functions
    # and the ratio of permeabilities.
    denominator_term = (sympy.cosh(a*d) + (mu_0 / mu) * sympy.sinh(a*d))**2

    # The force is in the x-direction (i_x).
    # We will represent this by creating a symbolic expression for the x-component.
    force_x_component = numerator_term / denominator_term

    # Print the final result in a clear, formatted way.
    print("The final expression for the force per unit area on the interface at x = d is:")
    
    # Use sympy's pretty printing for a more readable mathematical output.
    sympy.init_printing(use_unicode=True)
    
    # We create a "symbol" for the unit vector to make the output clear.
    i_x = sympy.Symbol("i_x")

    final_force_vector = force_x_component * i_x
    
    print("\nResult:")
    sympy.pprint(final_force_vector)

if __name__ == '__main__':
    display_force_equation()
