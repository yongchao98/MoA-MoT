def print_magnetic_field_expression():
    """
    This function prints the derived expression for the magnetic field H(x, z)
    for a stack of superconducting strips.

    The expression is valid under the condition |x| >> a, where the field from
    each strip is approximated by a magnetic dipole. The total field is the
    sum of the applied field and the field induced by the currents in the
    infinite stack of strips.
    """

    # The final derived expression for the magnetic field H(x, z) is printed below.
    # It includes the applied field 'H_a' and the induced field term.
    # The numbers '2' in the formula are explicitly included as requested.
    
    expression = "H(x, z) = H_a - ( (pi**2 * H0 * (w**2 - a**2)) / (2 * D**2) ) * ( (sinh(pi*x/D)**2 - sin(pi*z/D)**2) / (sinh(pi*x/D)**2 + sin(pi*z/D)**2)**2 )"
    
    print("The expression for the magnetic field H(x, z) is:")
    print(expression)

# Execute the function to print the final expression.
print_magnetic_field_expression()