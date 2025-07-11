def print_ac_loss_formula():
    """
    Prints the derived formula for the normalized AC loss in a superconductor.

    The problem asks for the loss per cycle per unit length (Q) for a superconductor
    in the critical state, with an elliptic cross-section, carrying a transport AC
    current with normalized amplitude i = Im/Ic, for i < 1.

    The result is presented in the standard normalized form: 2*pi*Q/(mu_0*Ic^2).
    """

    print("The normalized AC loss per cycle for a transport current in a superconductor with i < 1 is given by the following formula.")
    print("This formula is universal for any singly-connected conductor shape, including an ellipse of any aspect ratio.")
    print("-" * 80)

    # The formula is derived from the Norris model for transport current loss.
    # The final equation is: 2*pi*Q/(mu_0*Ic^2) = 2 * [ (1-i)ln(1-i) + (1+i)ln(1+i) - i^2 ]
    
    # We will print the equation components clearly.
    equation_lhs = "2 * pi * Q / (mu_0 * Ic^2)"
    # Using 'i**2' for i squared
    equation_rhs = "2 * [ (1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i**2 ]"

    print(f"Final Equation: {equation_lhs} = {equation_rhs}")
    
    print("-" * 80)
    print("where:")
    print("  Q:    Loss per cycle per unit length")
    print("  i:    Normalized current amplitude (Im/Ic)")
    print("  mu_0: Permeability of free space")
    print("  Ic:   Critical current of the bar")
    print("  ln:   Natural logarithm")
    print("  **:    Exponentiation (e.g., i**2 is i-squared)")

if __name__ == "__main__":
    print_ac_loss_formula()