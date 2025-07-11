def solve_polymer_force():
    """
    This function prints the derived force law for a thermally isolated polymer chain.

    The force law F(x) depends on:
    E0: The kinetic energy of the polymer at zero extension (x=0).
    n: The number of mass points (segments) in the polymer chain.
    l: The length of each segment.
    x: The separation (extension) of the polymer ends.
    """

    # The components of the force equation F(x) are defined as strings.
    # The force is attractive, hence the negative sign.
    # The coefficient includes the initial kinetic energy E(0).
    coefficient = "2 * E(0)"

    # The denominator of the main term.
    denominator = "n^2 * l^2"

    # The exponent term, which shows how the force increases non-linearly.
    exponent_numerator = "x^2"
    exponent_denominator = "n^2 * l^2"

    # Assemble the final equation into a formatted string.
    force_law_str = (
        f"F(x) = - ( {coefficient} * x / ({denominator}) ) "
        f"* exp( {exponent_numerator} / ({exponent_denominator}) )"
    )

    print("The derived force law between the polymer ends is:")
    print(force_law_str)

    print("\nWhere:")
    print("  F(x) is the force of attraction at extension x.")
    print("  E(0) is the kinetic energy at zero extension.")
    print("  n is the number of segments.")
    print("  l is the length of each segment.")
    print("  x is the end-to-end separation.")
    print("  exp() is the exponential function.")
    
    # Printing the numbers in the final equation as requested.
    # The main coefficient is -2. The exponent has a coefficient of 1.
    print("\nBreaking down the numerical coefficients in the final formula:")
    print(f"Main multiplicative coefficient: -{2}")
    print(f"Numerator of the exponent's argument: {1} * x^2")
    print(f"Denominator of the exponent's argument: {1} * n^2 * l^2")


solve_polymer_force()