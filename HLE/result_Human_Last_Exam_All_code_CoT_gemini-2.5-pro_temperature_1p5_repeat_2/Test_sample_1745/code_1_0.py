import sympy as sp

def display_edl_potential_expression():
    """
    This function defines and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution, psi(y),
    in a parallel-plate microchannel under the specified conditions.
    """
    # Define the symbolic variables used in the equation
    y, k, H, z1, beta = sp.symbols('y k H z_1 beta')
    
    # Define psi as a function of y
    psi = sp.Function('psi')(y)

    # Construct the final derived expression for the potential distribution
    # psi(y) = z_1 * (1 + beta * k) * sinh(k*(H/2 - y)) / sinh(k*H)
    expression = z1 * (1 + beta * k) * (sp.sinh(k * (H/2 - y)) / sp.sinh(k * H))

    # Create the full equation, psi(y) = expression
    final_equation = sp.Eq(psi, expression)

    # Print the final result in a human-readable format
    print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
    
    # sp.pprint is used for a clear, structured output of the symbolic equation
    sp.pprint(final_equation, use_unicode=True)

if __name__ == "__main__":
    display_edl_potential_expression()