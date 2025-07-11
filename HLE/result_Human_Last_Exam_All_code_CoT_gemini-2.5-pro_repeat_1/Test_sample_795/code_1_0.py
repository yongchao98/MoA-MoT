def print_magnetization_formula():
    """
    This function prints the derived analytical expression for the initial
    magnetization curve M(H) of the superconducting bar.

    The formula is valid for an applied field H in the range 0 <= H <= Jc * a.
    """
    # Define the variables for clarity in the formula string
    M_H = "M(H)"  # Magnetization as a function of H
    H = "H"      # Applied magnetic field
    Jc = "Jc"    # Critical current density
    a = "a"      # Half-width of the superconducting bar

    # Construct the formula string
    # The derived formula is: M(H) = H^2 / (2 * Jc * a) - H
    formula = f"{M_H} = {H}**2 / (2 * {Jc} * {a}) - {H}"

    # Print the final result
    print("The analytical expression for the initial magnetization curve is:")
    print(formula)
    print("\nWhere:")
    print("M(H) is the magnetization.")
    print("H is the applied magnetic field.")
    print("Jc is the critical current density.")
    print("a is the half-width of the bar in the x-direction.")
    print("\nThis expression is valid for 0 <= H <= Jc * a.")

# Execute the function to display the result
print_magnetization_formula()