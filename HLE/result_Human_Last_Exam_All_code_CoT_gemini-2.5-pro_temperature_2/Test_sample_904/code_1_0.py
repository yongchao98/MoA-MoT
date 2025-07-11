import sympy as sp

def solve_fluid_equation_coefficients():
    """
    This function uses symbolic mathematics to define and display the coefficients
    A(r) and B(r) of the governing linear differential equation for the fluid interface.

    The derivation is based on a force balance for a cylindrical system where the
    surface tension gamma is a function of the radial position r.
    """
    # Define symbols for the radial coordinate and functions for surface tension
    r = sp.Symbol('r')
    gamma = sp.Function('gamma')(r)

    # From the physical derivation, we identified the coefficients A(r) and B(r)
    # in the equation A(r)*xi'' + B(r)*xi' + C = 0.

    # Coefficient of the second derivative term (d^2(xi)/dr^2)
    A_r = gamma

    # Coefficient of the first derivative term (d(xi)/dr)
    # This includes a term from the non-constant surface tension
    B_r = gamma/r + sp.diff(gamma, r)

    # Print the equation form and the identified coefficients
    print("The governing linear equation has the form:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
    print("\nBased on the physical derivation, the coefficients are:")
    
    # "Output each number in the final equation" is interpreted as
    # printing the symbolic form of each coefficient term.
    print("\nCoefficient A(r):")
    sp.pprint(A_r)
    
    print("\nCoefficient B(r):")
    sp.pprint(B_r)

if __name__ == '__main__':
    solve_fluid_equation_coefficients()
