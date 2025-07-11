import sympy as sp

def solve_coefficients():
    """
    This function derives and prints the coefficients A(r) and B(r) for the
    governing linear equation of the fluid interface shape.
    """
    # Define the symbols used in the derivation.
    # r represents the radial position.
    # gamma represents the surface tension coefficient.
    r, gamma = sp.symbols('r gamma', positive=True)

    # From the linearized Young-Laplace equation in cylindrical coordinates,
    # the governing equation is:
    # gamma * d^2(xi)/dr^2 + (gamma/r) * d(xi)/dr - Delta_P_el(r) = 0

    # This is compared to the general form:
    # A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0

    # Coefficient of the second derivative term d^2(xi)/dr^2
    A_r = gamma

    # Coefficient of the first derivative term d(xi)/dr
    B_r = gamma / r

    # Print the derived coefficients
    print("For the governing linear equation of the interfacial shape xi(r):")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0\n")
    print("The coefficients A(r) and B(r) are:")
    sp.pprint(sp.Eq(sp.Function('A')(r), A_r), use_unicode=False)
    sp.pprint(sp.Eq(sp.Function('B')(r), B_r), use_unicode=False)

if __name__ == "__main__":
    solve_coefficients()