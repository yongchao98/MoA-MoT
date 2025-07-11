import math

def print_magnetic_field_expression():
    """
    This function prints the derived expression for the magnetic field H_z(x, z)
    for a stack of superconducting strips.
    """

    # The problem asks for the expression for the magnetic field H for Ha > H0 and |x| >> a.
    # The total field H_z is the sum of the applied field Ha and the field generated
    # by the currents in the superconducting stack, H_stack_z.
    #
    # H_z(x, z) = Ha + H_stack_z(x, z)
    #
    # Based on the plan outlined, the derived field contribution from the stack is:
    # H_stack_z = (pi * Jc * d * w**2 / (2 * D**2)) *
    #             (sinh(pi*x/D)**2 - sin(pi*z/D)**2) / (sinh(pi*x/D)**2 + sin(pi*z/D)**2)**2
    #
    # The following code prints this final expression term by term.

    print("The expression for the magnetic field component H_z(x, z) is:")
    print("\nH_z(x, z) = Ha + H_stack_z(x, z)\n")
    print("Where Ha is the applied magnetic field and H_stack_z is the field from the stack.\n")

    print("The expression for H_stack_z is:")
    # Printing the expression piece by piece to satisfy the "output each number" requirement
    # and to make it readable.
    print("H_stack_z(x, z) = ( (pi * Jc * d * w**2) / (2 * D**2) ) * NUMERATOR / DENOMINATOR**2\n")

    print("With the terms defined as:")
    print("pi = The mathematical constant pi (~3.14159)")
    print("Jc = Critical current density of the superconductor")
    print("d = Thickness of each strip")
    print("w = Half-width of each strip")
    print("D = Stacking interval between strips")
    print("x, z = Coordinates at which the field is calculated\n")

    print("And where:")
    print("NUMERATOR = sinh( (pi * x) / D )**2 - sin( (pi * z) / D )**2")
    print("DENOMINATOR = sinh( (pi * x) / D )**2 + sin( (pi * z) / D )**2\n")

    print("Combining all the terms, the full final expression is:")
    print("H_z(x, z) = Ha + ( (pi * Jc * d * w**2) / (2 * D**2) ) * (sinh((pi*x)/D)**2 - sin((pi*z)/D)**2) / (sinh((pi*x)/D)**2 + sin((pi*z)/D)**2)**2")

# Execute the function to print the result.
print_magnetic_field_expression()