def print_polymer_force_law():
    """
    Prints the derived force law for a thermally isolated, freely jointed polymer chain.

    The formula expresses the restoring force F(x) as a function of the
    end-to-end separation 'x', the number of segments 'n', the segment
    length 'l', and the polymer's kinetic energy at zero extension 'E(0)'.
    """

    # The derived mathematical formula for the force law.
    # The numbers 2 are explicitly included in the final printed equation as requested.
    force_law = "F(x) = (2 * E(0) * x / (n**2 * l**2)) * exp(-x**2 / (n**2 * l**2))"

    print("The force law for the thermally isolated polymer is:")
    print(force_law)

if __name__ == "__main__":
    print_polymer_force_law()