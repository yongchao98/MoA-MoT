def print_solution():
    """
    This function prints the closed-form expression for the given infinite product.
    The formula is derived using the Weierstrass product representation of the Gamma function.
    """

    print("The symbolic value of the infinite product is:")
    print("")
    print("  product_{n=3 to infinity} (1 - z^3 / n^3) = ")
    print("")
    print("                               1")
    print("  -------------------------------------------------------------------------")
    print("  (1 - z^3 / 1) * (1 - z^3 / 8) * Gamma(1-z) * Gamma(1-z*w) * Gamma(1-z*w^2)")
    print("")
    print("where 'w' is a complex cube root of unity, w = exp(i * 2 * pi / 3).")
    print("Gamma(x) is the Gamma function.")
    print("")
    print("The numbers present in the final equation are: 3, 3, 1, 3, 8, 1, 1, 2, 3.")

print_solution()