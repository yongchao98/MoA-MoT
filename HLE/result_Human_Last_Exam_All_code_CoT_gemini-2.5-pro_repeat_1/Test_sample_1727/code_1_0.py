import math

def solve_for_chi():
    """
    This function calculates the stable amplitude factor chi for the given geodynamic problem.
    The derivation shows that chi is a function of the dimensionless parameter H = k*h.
    Since H is not specified in the problem, we calculate the value for the canonical
    case H = 1, where the horizontal wavelength of the forcing is comparable to the
    layer thickness.

    The derived formula for chi is:
    chi(H) = (2 * H * cosh(H)) / (1 + (1 - 2 * H**2) * cosh(2*H))

    We will now compute this for H=1.
    """
    H = 1.0

    # Calculate the values of the hyperbolic functions
    cosh_H = math.cosh(H)
    cosh_2H = math.cosh(2 * H)

    # Calculate the numerator and the denominator of the expression for chi
    numerator = 2 * H * cosh_H
    denominator = 1 + (1 - 2 * H**2) * cosh_2H

    # Calculate the final value of chi
    chi = numerator / denominator

    # As requested, we will output the numbers that form the final equation for chi.
    print("The stable amplitude factor chi for H=1 is calculated as follows:")
    print(f"chi = (2 * {H} * cosh({H})) / (1 + (1 - 2 * {H}**2) * cosh(2 * {H}))")
    print(f"chi = (2 * {H} * {cosh_H:.5f}) / (1 + (1 - 2 * {H**2}) * {cosh_2H:.5f})")
    print(f"chi = ({numerator:.5f}) / (1 + ({-1.0}) * {cosh_2H:.5f})")
    print(f"chi = {numerator:.5f} / {denominator:.5f}")
    print(f"The final calculated value for chi is: {chi:.5f}")

solve_for_chi()