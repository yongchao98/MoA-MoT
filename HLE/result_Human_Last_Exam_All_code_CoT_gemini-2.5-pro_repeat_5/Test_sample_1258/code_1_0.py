import math

def print_demagnetizing_factor_expression():
    """
    Prints the analytical expression for the fluxmetric demagnetizing factor
    for a cylinder with uniform susceptibility chi=0.

    The expression is given in terms of a modulus 'k', where k^2 = 1 / (1 + g^2/4)
    and g is the length-to-diameter ratio of the cylinder.
    F(x) and E(x) are the complete elliptic integrals of the first and second kind,
    respectively, with modulus x.
    """

    expression_parts = [
        "N_f = (2 * k**2) / (pi * (1 - k**2) * sqrt(1 - k**2)) * ",
        "[( (1 + k**2) / 2 ) * F(sqrt(1 - k**2)) - E(sqrt(1 - k**2))]"
    ]

    print("The analytical expression for the fluxmetric demagnetizing factor (N_f) is:")
    print(" ".join(expression_parts))
    print("\nWhere:")
    print("  g = length-to-diameter ratio of the cylinder")
    print("  k^2 = 1 / (1 + g^2 / 4)")
    print("  F(x) is the complete elliptic integral of the first kind with modulus x.")
    print("  E(x) is the complete elliptic integral of the second kind with modulus x.")
    print("  pi is the mathematical constant {}...".format(math.pi))
    print("  sqrt(x) denotes the square root of x.")

if __name__ == '__main__':
    print_demagnetizing_factor_expression()
