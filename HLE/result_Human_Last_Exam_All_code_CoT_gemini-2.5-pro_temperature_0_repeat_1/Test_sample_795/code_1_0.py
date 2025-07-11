def print_magnetization_formula():
    """
    This function prints the derived analytical expression for the initial
    magnetization curve of a thin superconducting strip in a perpendicular
    magnetic field, based on the Bean critical-state model.
    """

    # The analytical formula for Magnetization (M) as a function of applied field (H)
    formula = "M(H) = H**2 / (2 * a * Jc) - H"

    # The range of validity for the formula
    validity_range = "0 <= H <= a * Jc"

    print("The analytical expression for the initial magnetization curve is:")
    print(f"  {formula}")
    print("\nThis expression is valid for the following range of the applied field H:")
    print(f"  {validity_range}")
    print("\nWhere the variables represent:")
    print("  M(H): Magnetization as a function of the applied field")
    print("  H: Applied magnetic field")
    print("  a: Half-width of the superconducting bar")
    print("  Jc: Critical current density (a constant)")
    print("  2: The number two in the denominator of the first term")

if __name__ == '__main__':
    print_magnetization_formula()