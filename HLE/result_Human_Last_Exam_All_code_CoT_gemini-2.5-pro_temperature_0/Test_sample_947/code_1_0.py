def print_magnetic_field_expression():
    """
    This function prints the derived expression for the magnetic field H_z(x, z)
    for a stack of superconducting strips under an applied field Ha.

    The expression is valid for |x| >> w, D, where w is the strip half-width
    and D is the stacking period.
    """

    # The expression for the z-component of the magnetic field
    hz_expression = "H_z(x, z) = Ha - (2*Jc*d/pi) * [sinh^2(pi*w/D) - sinh^2(pi*a/D)] * exp(-2*pi*|x|/D) * sin(2*pi*z/D)"

    # The expression relating the flux front 'a' to the applied field 'Ha'
    a_expression = "a = (D/pi) * arctanh[tanh(pi*w/D) * tanh(pi*Ha/(Jc*d))]"

    print("The expression for the magnetic field H_z for |x| >> a is given by:")
    print(hz_expression)
    print("\nIn this expression, the variables are:")
    print("  H_z(x, z): The magnetic field component in the z-direction at position (x, z).")
    print("  Ha: The applied external magnetic field.")
    print("  Jc: The critical current density of the superconductor.")
    print("  d: The thickness of each strip.")
    print("  w: The half-width of each strip (total width is 2w).")
    print("  D: The stacking period (distance between the centers of adjacent strips).")
    print("  pi: The mathematical constant pi (approx. 3.14159).")
    print("  exp, sin, sinh, arctanh, tanh: Standard mathematical functions.")
    print("  |x|: The absolute value of the x-coordinate.")
    print("  a: The position of the flux front, which depends on Ha.")
    print("\nThe flux front position 'a' is determined by the applied field 'Ha' through the following relation:")
    print(a_expression)

if __name__ == '__main__':
    print_magnetic_field_expression()