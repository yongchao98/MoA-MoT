def print_potential_distribution_formula():
    """
    This function constructs and prints the final expression for the
    Electrical Double-Layer (EDL) potential distribution (Psi(y)).
    """
    
    # Define the components of the formula as string variables
    psi_y = "Ψ(y)"
    
    # Bottom surface zeta potential (z1)
    zeta_1 = "z_1"
    
    # Slip-dependant part with slip length (beta) and Debye-Huckel parameter (k)
    slip_component = "(1 + β*k)"
    
    # Hyperbolic sine part of the formula
    # H is the channel height, y is the vertical position
    sinh_part = "sinh(k*(H - y)) / sinh(k*H)"
    
    # Construct the final formula string
    # The equation describes the potential Psi at a vertical position y.
    # z_1: Zeta potential of the bottom surface.
    # beta: Slip length.
    # k: Debye-Huckel parameter.
    # H: Height of the microchannel.
    final_formula = f"{psi_y} = {zeta_1} * {slip_component} * [ {sinh_part} ]"

    # Print the final formatted equation
    print("The expression for the Electrical double-layer potential distribution Ψ(y) is:")
    print(final_formula)

if __name__ == "__main__":
    print_potential_distribution_formula()