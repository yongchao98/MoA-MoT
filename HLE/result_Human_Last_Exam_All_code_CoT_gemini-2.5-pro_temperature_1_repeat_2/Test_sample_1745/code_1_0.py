def solve_potential_distribution():
    """
    This function prints the derived expression for the electrical double-layer (EDL)
    potential distribution psi(y) in the microchannel.
    """

    # The problem asks for the mathematical expression as the output.
    # The derivation leads to the following formula for the potential psi(y),
    # where z1 is the zeta potential of the bottom surface,
    # beta is the slip length, k is the Debye-Huckel parameter,
    # H is the channel height, and y is the vertical coordinate.

    expression = (
        "The expression for the Electrical double-layer potential distribution psi(y) is:\n\n"
        "         z1 * (1 + beta * k) * sinh(k * (H/2 - y))\n"
        "psi(y) = ---------------------------------------------\n"
        "                       sinh(k * H)                     "
    )

    print(expression)

    # To explicitly show the components of the equation as requested:
    print("\nWhere the terms in the final equation are:")
    print("psi(y): Electrical double-layer potential at position y")
    print("z1: Zeta potential of the bottom surface (j=1)")
    print("beta: Slip length")
    print("k: Debye-Huckel parameter")
    print("H: Height of the microchannel")
    print("y: Coordinate perpendicular to the plates (-H/2 <= y <= H/2)")
    print("sinh: Hyperbolic sine function")

if __name__ == "__main__":
    solve_potential_distribution()