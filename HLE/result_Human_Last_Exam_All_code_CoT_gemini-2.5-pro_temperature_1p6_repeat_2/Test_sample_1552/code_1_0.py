from fractions import Fraction

def find_second_heat_kernel_coefficient():
    """
    This function calculates and displays the second Seeley-DeWitt coefficient (a_2)
    for a massless gauged Dirac spinor field in 4 dimensions.

    The coefficient is an expression involving geometric and gauge field invariants.
    N represents the dimension of the gauge group representation.
    """

    # The coefficients are derived from established results in the literature for the
    # traced heat kernel density, tr(a_2(D^2)). The trace is over both spinor
    # (dimension 4) and gauge group (dimension N) indices.
    
    # Coefficient for the square of the Riemann tensor, R_abcd * R^{abcd}
    c_riemann = Fraction(4, 1) * Fraction(-7, 360)

    # Coefficient for the square of the Ricci tensor, R_ab * R^{ab}
    c_ricci = Fraction(4, 1) * Fraction(-1, 45)

    # Coefficient for the square of the Ricci scalar, R^2
    c_scalar = Fraction(4, 1) * Fraction(1, 72)

    # Coefficient for the square of the gauge field strength, Tr_G(F_ab * F^{ab})
    # The -1/2 comes from the standard formula, the 4 is from the trace over spinor indices.
    c_gauge = Fraction(-1, 2) * Fraction(4, 1)

    print("The second heat kernel coefficient density, traced over the fiber, is given by:")
    print("a_2(x) = (C1 * N) * R_{\\mu\\nu\\rho\\sigma}R^{\\mu\\nu\\rho\\sigma} + (C2 * N) * R_{\\mu\\nu}R^{\\mu\\nu} + (C3 * N) * R^2 + C4 * Tr_G(F_{\\mu\\nu}F^{\\mu\\nu})\n")
    print("Where N is the dimension of the gauge representation, and the coefficients are:")
    
    print(f"C1 = {c_riemann.limit_denominator()} = {c_riemann.numerator}/{c_riemann.denominator}")
    print(f"C2 = {c_ricci.limit_denominator()} = {c_ricci.numerator}/{c_ricci.denominator}")
    print(f"C3 = {c_scalar.limit_denominator()} = {c_scalar.numerator}/{c_scalar.denominator}")
    print(f"C4 = {c_gauge.limit_denominator()} = {c_gauge.numerator}/{c_gauge.denominator}\n")

    print("The final expression is:")
    print(f"a_2(x) = ({c_riemann.numerator}/{c_riemann.denominator}*N) * R_{\\mu\\nu\\rho\\sigma}R^{\\mu\\nu\\rho\\sigma} + ({c_ricci.numerator}/{c_ricci.denominator}*N) * R_{\\mu\\nu}R^{\\mu\\nu} + ({c_scalar.numerator}/{c_scalar.denominator}*N) * R^2 + ({c_gauge.numerator}/{c_gauge.denominator}) * Tr_G(F_{\\mu\\nu}F^{\\mu\\nu})")


find_second_heat_kernel_coefficient()