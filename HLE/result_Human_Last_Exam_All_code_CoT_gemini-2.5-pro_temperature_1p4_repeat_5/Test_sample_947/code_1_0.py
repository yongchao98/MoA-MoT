def print_magnetic_field_expression():
    """
    This function prints the derived expressions for the magnetic field components H_x and H_z.
    """

    # Define the expressions for the magnetic field components as formatted strings.
    # Note: sinh is the hyperbolic sine, and sin is the standard sine.
    # H0, w, D, Ha are parameters from the problem description.
    # pi is the mathematical constant.
    # x and z are the coordinates where the field is calculated.

    H_x_expression = (
        "H_x(x, z) = (H0 * pi**2 * w**2 / (4 * D**2)) * "
        "(sin(2*pi*z/D) * sinh(2*pi*x/D)) / (sin(pi*z/D)**2 + sinh(pi*x/D)**2)**2"
    )

    H_z_expression = (
        "H_z(x, z) = Ha - (H0 * pi**2 * w**2 / (2 * D**2)) * "
        "(sin(pi*z/D)**2 - sinh(pi*x/D)**2) / (sin(pi*z/D)**2 + sinh(pi*x/D)**2)**2"
    )

    # Print the final expressions for the magnetic field vector H = (H_x, 0, H_z).
    print("The expression for the magnetic field H = (H_x, 0, H_z) is given by its components:")
    print("\nComponent H_x:")
    print(H_x_expression)
    print("\nComponent H_z:")
    print(H_z_expression)

# Execute the function to print the result.
print_magnetic_field_expression()