def find_computational_factor():
    """
    This function identifies and prints the value of the computational factor (mushy zone constant)
    from the prior, simulation-only work that established the enthalpy-porosity method for melting.

    The value is taken from the 1988 paper by A. D. Brent, V. R. Voller, and K. J. Reid,
    "A numerical study of the melting of a pure metal from an isothermal wall".
    """

    # The computational factor, often denoted as 'C', is the mushy zone constant.
    # In the original 1988 paper, this value was set to 1.6 x 10^6.
    coefficient = 1.6
    base = 10
    exponent = 6

    # Calculate the final value
    value = coefficient * (base ** exponent)

    # Print the result in the format of the final equation
    print("The computational factor 'C' from the prior simulation-only work is given by:")
    print(f"C = {coefficient} * {base}^{exponent}")
    print(f"C = {int(value):,}") # Print with commas for readability

find_computational_factor()