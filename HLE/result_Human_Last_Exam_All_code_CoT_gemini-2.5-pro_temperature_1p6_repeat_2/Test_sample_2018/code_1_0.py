def solve_historical_value():
    """
    This function identifies and prints the value of the computational factor
    from the original paper describing the enthalpy-porosity method for
    phase-change simulation.
    """

    # The value is given in scientific notation.
    # The value comes from the 1987 paper by Voller and Prakash.
    mantissa = 1.6
    base = 10
    exponent = 6

    # Calculate the full value
    value = mantissa * (base ** exponent)

    # Print the equation as requested
    print("The computational factor (C) from the prior work is:")
    print(f"{mantissa} * {base}^{exponent} = {value:e}")

solve_historical_value()