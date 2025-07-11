import math

def display_ac_loss_formula():
    """
    This function prints the normalized AC loss formula for a superconductor
    with an elliptical cross-section carrying a transport current.
    """

    # The formula is derived from the critical-state model (Norris's solution).
    
    # Define the components of the formula as strings for clear printing.
    lhs = "2 * pi * Q / (mu_0 * Ic^2)"
    rhs = "2 * ((1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i^2)"
    
    print("The normalized AC loss per cycle per unit length is given by the formula:")
    print(f"\n{lhs} = {rhs}\n")
    
    print("Where:")
    print("  Q: Loss per cycle per unit length.")
    print("  mu_0: Vacuum permeability.")
    print("  Ic: Critical current of the bar.")
    print("  i: Normalized current amplitude, i = Im / Ic, where Im is the current amplitude.")
    print("  ln: Natural logarithm.")
    print("\nThis formula is valid for i < 1 and holds for any elliptical cross-section.")

# Execute the function to print the formula.
display_ac_loss_formula()