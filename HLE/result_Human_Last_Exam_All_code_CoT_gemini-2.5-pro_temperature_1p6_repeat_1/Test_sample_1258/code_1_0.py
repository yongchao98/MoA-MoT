def print_demagnetizing_factor_expression():
    """
    This function prints the analytical expression for the fluxmetric
    demagnetizing factor of a cylinder for the case of zero susceptibility.

    The expression is based on the magnetometric demagnetizing factor for a
    uniformly magnetized cylinder, as these two factors are equal when chi=0.

    Variables used in the expression:
    N: The demagnetizing factor (dimensionless).
    g: The length-to-diameter ratio of the cylinder (L/D).
    k: The modulus parameter, defined as k^2 = 1 / (1 + g^2 / 4).
    pi: The mathematical constant pi.
    F(x): The complete elliptic integral of the first kind with modulus x.
    E(x): The complete elliptic integral of the second kind with modulus x.
    """

    expression = (
        "The analytical expression for the fluxmetric demagnetizing factor N is:\n\n"
        "N = 2/g**2 - 8/(pi * g**3 * k) * [F(sqrt(1-k**2)) - E(sqrt(1-k**2))]"
    )

    print(expression)

if __name__ == "__main__":
    print_demagnetizing_factor_expression()
