def print_magnetic_field_formula():
    """
    This function prints the derived formula for the magnetic field of a stack
    of superconducting strips.
    """

    # The formula is derived based on the physical model described.
    # It represents the z-component of the magnetic field H_z at a point (x, z).

    print("The final expression for the z-component of the magnetic field H_z(x, z) is given by:")
    print("")
    print("H_z(x, z) = Ha + H_induced_z(x, z)")
    print("")
    print("where Ha is the applied magnetic field, and the induced field is:")
    print("")
    # Using the result derived from Brandt's formulation for an array of magnetic dipoles.
    # H_induced_z = (m_z / (2*D**2)) * (sin(pi*z/D)**2 - sinh(pi*x/D)**2) / (sinh(pi*x/D)**2 + sin(pi*z/D)**2)**2
    # Substitute m_z = -H0 * pi * w**2 and H0 = Jc * d / pi --> m_z = -Jc * d * w**2
    
    numerator_factor = "Jc * d * w**2"
    denominator_factor = "2 * D**2"
    
    numerator_terms = "sinh(pi*x/D)**2 - sin(pi*z/D)**2"
    denominator_terms = "(sinh(pi*x/D)**2 + sin(pi*z/D)**2)**2"

    print("                " + numerator_factor + "      " + numerator_terms)
    print("H_induced_z = -  -------------------- * -------------------------------------------")
    print("                   " + denominator_factor + "          " + denominator_terms)
    print("")
    print("Full equation:")
    print("                                      Jc * d * w**2      sinh(pi*x/D)**2 - sin(pi*z/D)**2")
    print("H_z(x, z) = Ha  -  ( ----------------- ) * ( ----------------------------------------- )")
    print("                                         2 * D**2          (sinh(pi*x/D)**2 + sin(pi*z/D)**2)**2")
    print("")
    print("Variables:")
    print("  Ha: Applied magnetic field in the z-direction.")
    print("  Jc: Critical current density of the superconductor.")
    print("  d: Thickness of each strip.")
    print("  w: Half-width of each strip (total width is 2w).")
    print("  D: Stacking interval between the strips.")
    print("  x, z: Coordinates where the field is calculated.")
    print("  pi: The mathematical constant pi (approx. 3.14159).")

print_magnetic_field_formula()