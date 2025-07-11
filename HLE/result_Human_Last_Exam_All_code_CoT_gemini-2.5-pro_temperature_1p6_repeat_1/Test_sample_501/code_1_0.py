def display_force_law():
    """
    This function constructs and prints the derived force law for a thermally
    isolated freely jointed polymer chain.
    """

    # The derived force law is F(x) = (2 * E(0) * x) / (n^2 * l^2)
    # We will define the numerical constants present in this equation.
    coefficient = 2
    n_exponent = 2
    l_exponent = 2

    # The variables in the equation are:
    # F(x): The force of attraction between the ends at extension x.
    # E(0): The kinetic energy of the polymer at zero extension (x=0).
    # x: The separation of the polymer ends.
    # n: The number of mass points (segments).
    # l: The length of each segment.

    # Now, we print the final equation in a clear, readable format.
    print("The derived force law between the polymer ends for a thermally isolated chain is:")
    print(f"F(x) = ({coefficient} * E(0) * x) / (n^{n_exponent} * l^{l_exponent})")
    print("\nThis result is valid for small extensions x.")

# Execute the function to display the result.
display_force_law()