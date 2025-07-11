def print_ac_loss_formula():
    """
    This function prints the normalized AC loss formula for a superconducting
    elliptic bar carrying a transport current, valid for i = Im/Ic < 1.
    """

    # The expression is derived from the analytical solution by W. T. Norris (1970).
    # Norris's formula for the loss per cycle per unit length, Q, is:
    # Q = (mu_0 * Ic**2 / pi) * F(i)
    # where F(i) = [ (1 - i/2)*ln(1-i) + (1 + i/2)*ln(1+i) - i**2 / 2 ]

    # The requested normalized form is 2*pi*Q / (mu_0 * Ic**2).
    #
    # Performing the normalization:
    # (2 * pi * Q) / (mu_0 * Ic**2) = (2 * pi / (mu_0 * Ic**2)) * [ (mu_0 * Ic**2 / pi) * F(i) ]
    #                              = 2 * F(i)
    
    # Define the Left-Hand Side (LHS) and Right-Hand Side (RHS) of the equation.
    normalized_loss_lhs = "2*pi*Q/(μ_0 * Ic^2)"
    normalized_loss_rhs = "2 * [ (1 - i/2)*ln(1 - i) + (1 + i/2)*ln(1 + i) - i**2/2 ]"

    print("The formula for the normalized AC loss in an elliptic superconductor is:")
    print("-" * 60)
    print(f"{normalized_loss_lhs} = {normalized_loss_rhs}")
    print("-" * 60)
    print("\nWhere:")
    print("  Q    = Loss per cycle per unit length")
    print("  μ_0  = Vacuum permeability")
    print("  Ic   = Critical current of the bar")
    print("  i    = Im / Ic, the normalized current amplitude (for i < 1)")
    print("  Im   = Amplitude of the transport AC current")
    print("  ln() = Natural logarithm")
    print("\nThis formula is valid for an elliptic bar with any aspect ratio (a/b).")

# Execute the function to print the formula.
print_ac_loss_formula()