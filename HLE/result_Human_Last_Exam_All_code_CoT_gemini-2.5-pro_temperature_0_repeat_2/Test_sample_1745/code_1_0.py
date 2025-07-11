def print_potential_distribution():
    """
    Prints the derived expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y).
    """

    # Define the variables in the equation for clarity
    # psi(y): EDL potential at position y
    # z_1: Zeta potential of the bottom surface (j=1)
    # beta: Slip length
    # k: Debye-Huckel parameter
    # H: Height of the microchannel
    # y: Position variable across the channel height (-H/2 to H/2)
    # sinh: Hyperbolic sine function

    expression = "psi(y) = z_1 * (1 + beta * k) * sinh(k * (H/2 - y)) / sinh(k * H)"

    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    print(expression)

if __name__ == "__main__":
    print_potential_distribution()