def print_capacitance_formula():
    """
    This script prints the derived formula for the gate capacitance per unit area (C_g)
    for the given quantum Hall system.
    It explicitly shows the numerical coefficients in the final equation as requested.
    """

    # The derived formula is C_g = (4 * e^2 * B) / (h * V_1).
    # The numerical values in this formula are:
    numerator_coefficient = 4
    power_of_e = 2
    h_coefficient = 1
    v1_coefficient = 1

    # Printing the final equation with all numbers explicitly shown.
    # Here, V_1 represents the variable V1 from the problem.
    print("The formula for the gate capacitance per unit area (C_g) is:")
    print(f"C_g = ({numerator_coefficient} * e^{power_of_e} * B) / ({h_coefficient} * h * {v1_coefficient} * V_1)")

# Execute the function to print the result.
print_capacitance_formula()