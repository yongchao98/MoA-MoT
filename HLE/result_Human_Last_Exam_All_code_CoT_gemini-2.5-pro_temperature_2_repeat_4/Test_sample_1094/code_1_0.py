import sympy

def solve_ac_loss():
    """
    This function derives and prints the normalized AC loss for an elliptical
    superconductor based on the Norris model.
    """
    # Define symbols
    # i represents the normalized current Im/Ic
    i = sympy.Symbol('i')
    
    # Norris formula for AC loss per cycle per unit length (Q) for i < 1
    # We represent the core part of the formula here.
    # The full formula is Q = (mu_0 * Ic**2 / pi) * loss_function
    loss_function = (1 - i) * sympy.ln(1 - i) + (1 + i) * sympy.ln(1 + i) - i**2
    
    # The user asks for the normalized loss: 2*pi*Q / (mu_0 * Ic**2)
    # Substituting Q: (2*pi / (mu_0 * Ic**2)) * (mu_0 * Ic**2 / pi) * loss_function
    # The prefactors simplify to just '2'.
    normalized_loss = 2 * loss_function
    
    # Format the final expression for printing
    final_expression = f"2 * [ (1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i**2 ]"

    print("The standard form for AC losses, 2*pi*Q/(mu_0*Ic^2), as a function of i = Im/Ic is:")
    print(final_expression)

# Execute the function to print the result
solve_ac_loss()