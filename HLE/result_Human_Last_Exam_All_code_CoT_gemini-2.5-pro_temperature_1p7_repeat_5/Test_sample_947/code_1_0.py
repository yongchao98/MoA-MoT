def solve_superconducting_stack_field():
    """
    This function calculates and prints the expression for the magnetic field
    of a stack of superconducting strips under an applied field Ha > H0.
    The expression is valid in the far-field limit (|x| >> w).
    """

    # Define the components of the final expression as strings
    field_H_z = "H_z(x, z)"
    applied_field = "Ha"

    # Define the magnetic moment 'm'
    # m = pi * H0 * w^2 * tanh^2(Ha/H0)
    # where H0 = Jc * d / pi
    moment_m_definition = "m = pi * H0 * w**2 * tanh**2(Ha/H0)"
    h0_definition = "H0 = Jc * d / pi"

    # Define the geometric factor from the dipole stack summation
    # F(x,z) = (pi / D**2) * [ sinh(2*pi*x/D) * sin(2*pi*z/D) / (cosh(2*pi*x/D) - cos(2*pi*z/D))**2 ]
    geometric_factor_numerator = "sinh(2*pi*x/D) * sin(2*pi*z/D)"
    geometric_factor_denominator = "(cosh(2*pi*x/D) - cos(2*pi*z/D))**2"
    induced_field_term = f"m * (pi / D**2) * ( {geometric_factor_numerator} / {geometric_factor_denominator} )"

    # Construct the final equation for H_z(x, z)
    final_equation = f"{field_H_z} = {applied_field} - {induced_field_term}"

    # Print the results in a structured way
    print("The expression for the z-component of the magnetic field H_z(x, z) is:")
    print(final_equation)
    print("\nwhere the magnetic moment per unit length 'm' is given by:")
    print(moment_m_definition)
    print("\nand the characteristic field 'H0' is:")
    print(h0_definition)

# Execute the function to print the solution
solve_superconducting_stack_field()