import sympy as sp

def display_coefficients():
    """
    Identifies and prints the coefficients A(r) and B(r) for the governing
    linear equation of a fluid interface under surface tension and electric field.
    """
    # Define symbols for the variables in the problem
    r = sp.Symbol('r')         # Radial position
    gamma = sp.Symbol('gamma') # Surface tension

    # The governing equation is given in the form:
    # A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0
    #
    # The terms involving derivatives of xi(r) arise from the linearized
    # Young-Laplace equation for surface tension pressure in cylindrical coordinates.
    # The surface tension pressure term is: gamma * (d^2(xi)/dr^2 + (1/r) * d(xi)/dr)
    #
    # By comparing the physical equation with the given general form,
    # we can identify the coefficients A(r) and B(r).

    # Coefficient of the second derivative, d^2(xi)/dr^2
    A_r = gamma

    # Coefficient of the first derivative, d(xi)/dr
    B_r = gamma / r

    print("For the governing equation:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0\n")

    print("The coefficient A(r) is the coefficient of the second derivative term.")
    print("A(r) =", A_r)
    print("This is the constant surface tension, gamma.\n")

    print("The coefficient B(r) is the coefficient of the first derivative term.")
    print("B(r) =", B_r)
    print("This is the surface tension, gamma, divided by the radial position, r.")

if __name__ == '__main__':
    display_coefficients()