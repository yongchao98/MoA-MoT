import math

def print_magnetic_field_expression():
    """
    This function prints the derived expression for the magnetic field
    components H_x and H_z for a stack of superconducting strips.
    """

    # The problem asks for the expression for the magnetic field H.
    # The magnetic field is a vector with components H_x and H_z (H_y = 0 by symmetry).
    # H(x, z) = H_x(x, z) * x_hat + H_z(x, z) * z_hat
    # The total field is H_total = H_applied + H_induced.
    # H_applied is (0, 0, Ha).

    # The common coefficient for the induced field, derived from the physics of the system.
    # It contains the magnetic moment of a single strip and geometric factors of the stack.
    # Note: we use "pi" to represent the mathematical constant π.
    # Note: sgn(x) is the sign function, which is -1 for x<0, 0 for x=0, and 1 for x>0.
    # Note: |x| is the absolute value of x.
    # The number '2' appears multiple times from the derivation.
    
    coefficient = "(2 * pi**2 * H0 * w**2 / D**2) * (1 - 2*H0/Ha)"

    # Expression for the x-component of the magnetic field
    h_x_expression = f"H_x(x, z) = sgn(x) * {coefficient} * exp(-2*pi*|x|/D) * sin(2*pi*z/D)"

    # Expression for the z-component of the magnetic field
    # It is the sum of the applied field Ha and the induced field.
    h_z_expression = f"H_z(x, z) = Ha - {coefficient} * exp(-2*pi*|x|/D) * cos(2*pi*z/D)"

    print("The derived expressions for the magnetic field components are:")
    print("-" * 60)
    print(h_x_expression)
    print("-" * 60)
    print(h_z_expression)
    print("-" * 60)
    print("Where:")
    print("  Ha: Applied magnetic field")
    print("  H0: Full penetration field for a single strip (H0 = Jc*d/pi)")
    print("  Jc: Critical current density")
    print("  w: Half-width of the strips (total width is 2w)")
    print("  d: Thickness of the strips")
    print("  D: Stacking interval")
    print("  (x, z): Coordinates of the point of interest")
    print("  pi: The mathematical constant " + str(math.pi) + "...")

# Execute the function to print the result
print_magnetic_field_expression()
