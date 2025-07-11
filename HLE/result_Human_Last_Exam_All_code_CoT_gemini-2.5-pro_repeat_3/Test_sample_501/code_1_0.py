def solve_polymer_force():
    """
    This function constructs and prints the final force law for a thermally
    isolated freely jointed polymer chain.

    The variables in the equation are:
    f(x): The force as a function of extension x.
    E(0): The kinetic energy of the polymer at zero extension (x=0).
    x:    The separation distance between the polymer ends.
    n:    The number of links in the polymer chain.
    l:    The length of each link (using 'l' for \ell).
    exp(): The exponential function.
    """

    # Construct the components of the equation as strings
    coefficient = "(2 * E(0) * x / (n**2 * l**2))"
    exponent = "(x**2 / (n**2 * l**2))"

    # Assemble the final equation string
    force_law = f"f(x) = -{coefficient} * exp{exponent}"

    # Print the final result
    print("The derived force law for the thermally isolated polymer is:")
    print(force_law)

solve_polymer_force()