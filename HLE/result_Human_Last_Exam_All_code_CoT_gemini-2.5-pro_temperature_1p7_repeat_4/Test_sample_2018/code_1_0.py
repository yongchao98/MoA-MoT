def solve_computational_factor():
    """
    This function determines and prints the value of the computational factor
    used in the prior published work for the enthalpy-porosity method.

    The value is identified from the 1987 paper by Voller and Prakash as 1.6 x 10^6.
    """

    # The computational factor is often called the mushy zone constant, C.
    # Base of the number in scientific notation.
    base_value = 1.6
    # Exponent of the number in scientific notation.
    exponent_value = 6

    # Calculate the full value.
    final_value = base_value * (10 ** exponent_value)

    # Print the equation as requested, showing each number.
    print(f"The value for the computational factor in the prior published work was:")
    print(f"{base_value} * 10^{exponent_value} = {final_value:n}")

solve_computational_factor()