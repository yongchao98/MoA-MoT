def print_magnetic_field_expression():
    """
    This function prints the mathematical expression for the z-component of the
    magnetic field Hz(x, z) for a stack of superconducting strips.
    """

    print("The expression for the z-component of the magnetic field, Hz(x, z), is given by:")
    print("--------------------------------------------------------------------------------")
    print("Hz(x, z) = Ha + H_induced_z\n")
    print("The induced field, H_induced_z, is found using complex analysis:")
    print("H_induced_z = (H0 / 2) * [atan2(Si, Sr - C) - atan2(Si, Sr)]\n")
    print("where 'atan2(y, x)' is the two-argument arctangent function.\n")
    print("The terms in the equation are defined as follows:\n")
    print("--- Input Parameters and Constants ---")
    print("Ha: Applied magnetic field in the z-direction.")
    print("H0: Full penetration field of a single strip, defined as Jc * d / pi.")
    print("Jc: Critical current density of the superconductor.")
    print("d: Thickness of each strip.")
    print("w: Half-width of each strip (strip extends from x = -w to x = +w).")
    print("D: Stacking interval (distance between the centers of adjacent strips).")
    print("pi: The mathematical constant pi (approx. 3.14159).\n")
    print("--- Intermediate Variables ---")
    print("The calculation depends on three intermediate variables: C, Sr, and Si.\n")
    print("1. C is a constant that depends on the geometry of the strips:")
    print("   C = sin(pi * w / D)^2\n")
    print("2. Sr is the real part of an intermediate complex variable:")
    print("   Sr = 0.5 * (1 - cos(2*pi*x/D) * cosh(2*pi*z/D))\n")
    print("3. Si is the imaginary part of the same intermediate complex variable:")
    print("   Si = 0.5 * sin(2*pi*x/D) * sinh(2*pi*z/D)\n")
    print("Here, 'cos' and 'sin' are the standard trigonometric functions, and")
    print("'cosh' and 'sinh' are the hyperbolic functions.")
    print("--------------------------------------------------------------------------------")

if __name__ == "__main__":
    print_magnetic_field_expression()
