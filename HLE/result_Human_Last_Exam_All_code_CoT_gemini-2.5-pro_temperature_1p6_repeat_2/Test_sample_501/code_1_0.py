def print_force_law():
    """
    This function prints the derived force law for a thermally isolated polymer chain.
    The force 'f' is expressed in terms of the initial kinetic energy 'E(0)',
    the number of segments 'n', the length of each segment 'l', and the
    end-to-end separation 'x'.
    """

    # Define symbolic variables as strings for clear output
    E_0 = "E(0)"  # Kinetic energy at zero extension
    n = "n"       # Number of segments
    l = "ell"     # Length of each segment
    x = "x"       # Separation of the ends

    # The derived force law for small x is f = - (2 * E(0) * x) / (n^2 * l^2)
    # The only numerical constant in the final equation is 2.
    numerator_constant = 2

    # Construct and print the final equation
    print("The force law f(x) for a thermally isolated polymer chain at small extension x is:")
    print(f"f(x) = - ( {numerator_constant} * {E_0} * {x} ) / ( {n}**2 * {l}**2 )")

# Execute the function to display the result
print_force_law()