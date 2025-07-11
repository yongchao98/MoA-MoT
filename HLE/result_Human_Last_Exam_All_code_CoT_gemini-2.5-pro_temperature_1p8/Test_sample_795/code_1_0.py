import math

def print_magnetization_formula():
    """
    This function prints the analytical expression for the initial magnetization curve
    of a superconducting slab in the Bean critical-state model.
    """

    # The equation involves variables representing physical quantities and one constant number.
    # Variable definitions:
    # M: Magnetization
    # H: Applied magnetic field strength
    # Jc: Critical current density
    # a: Half-width of the superconducting bar

    # The constant number in the equation is 2.
    constant_number = 2

    # Constructing the expression string.
    # The formula describes the magnetization M as a function of the applied field H.
    # The term H**2 represents H squared.
    expression = f"M(H) = H**2 / ({constant_number} * Jc * a) - H"
    
    # This expression is valid up to the full penetration field, H_p = Jc * a.
    validity_condition = "This expression is valid for applied fields in the range: 0 <= H <= Jc * a"

    # Print the final result in a clear and understandable format.
    print("The analytical expression for the initial magnetization curve is:")
    print(expression)
    print("\n" + "="*50)
    print("Where:")
    print("  M(H): Magnetization as a function of H")
    print("  H:    Magnitude of the applied magnetic field")
    print("  Jc:   The constant critical-current density")
    print("  a:    The half-width of the superconducting bar")
    print(f"  {constant_number}:    The number '{constant_number}' is a constant factor in the denominator.")
    print("="*50 + "\n")
    print(validity_condition)

# Execute the function to display the result.
print_magnetization_formula()
