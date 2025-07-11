def find_computational_factor():
    """
    This function identifies and prints the value of the computational factor
    from the original simulation-only paper on the enthalpy-porosity method.
    """

    # The computational factor is a constant C, representing the mushy zone.
    # It is expressed in scientific notation.
    mantissa = 1.6
    exponent = 6

    print("The computational factor, C, used in the Carman-Kozeny source term can be found in the original paper describing the method.")
    print("In the prior, simulation-only work by Voller and Prakash (1987), this value was set for the example calculations.")
    print("\nThe final equation for the factor C is:")
    
    # Per the instructions, output each number in the final equation.
    # The numbers in the equation C = 1.6 * 10^6 are 1.6 and 6.
    print(f"C = {mantissa} * 10^{exponent}")
    
    final_value = mantissa * (10**exponent)
    print(f"\nThis evaluates to: {int(final_value)}")

find_computational_factor()