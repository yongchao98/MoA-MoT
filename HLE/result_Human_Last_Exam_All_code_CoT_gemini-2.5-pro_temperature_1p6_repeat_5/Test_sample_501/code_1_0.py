def print_force_law():
    """
    Prints the derived force law for a thermally isolated,
    freely jointed polymer chain.
    """
    
    print("The force law F(x) between the ends of the thermally isolated polymer is given by the following equation:")
    # The equation includes the numeric constants 3 and 2 as derived.
    force_equation = "F(x) = (3 * E(0) * x) / (n^2 * l^2) * exp( - (3 * x^2) / (2 * n^2 * l^2) )"
    print("\n" + force_equation + "\n")
    print("Where:")
    print("  F(x) is the force of attraction between the polymer ends.")
    print("  x    is the separation of the ends.")
    print("  E(0) is the kinetic energy of the polymer at zero extension (x=0).")
    print("  n    is the number of segments in the polymer chain (assumed to be large).")
    print("  l    is the length of a single segment.")
    print("  exp() is the exponential function e^().")

if __name__ == '__main__':
    print_force_law()