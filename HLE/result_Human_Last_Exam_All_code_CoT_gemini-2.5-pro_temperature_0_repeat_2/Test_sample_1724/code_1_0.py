def print_frequency_correction_coefficients():
    """
    This function prints the numerical coefficients found in the derived formula for
    the second-order frequency correction, omega_2.

    The formula is:
    omega_2 = - (sqrt(A*gamma) / B) * (C*gamma^2 + D*gamma + E)
    """

    # Coefficients from the derived formula
    A = 3
    B = 16
    C = 6
    D = 5
    E = 14

    print("The derived expression for the frequency correction omega_2 is of the form:")
    print("omega_2 = - (sqrt(A*gamma) / B) * (C*gamma^2 + D*gamma + E)")
    print("\nThe numerical coefficients in this equation are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")

print_frequency_correction_coefficients()