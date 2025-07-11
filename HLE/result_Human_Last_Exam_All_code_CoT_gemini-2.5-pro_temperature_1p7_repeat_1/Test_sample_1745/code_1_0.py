import sympy as sp

def solve_edl_potential():
    """
    This function defines and prints the symbolic expression for the 
    Electrical Double-Layer (EDL) potential distribution psi(y)
    for the given boundary conditions.
    """
    # Define the symbols used in the equation
    y, z_1, beta, k, H = sp.symbols('y z_1 beta k H')
    psi = sp.Function('psi')(y)

    # Define the slip-dependent zeta potential at the bottom surface (j=1)
    # The problem states z_2 = 0, so zeta_2 is 0.
    zeta_1 = z_1 * (1 + beta * k)

    # Construct the final derived expression for psi(y)
    # psi(y) = zeta_1 * sinh(k*(H/2 - y)) / sinh(k*H)
    numerator = sp.sinh(k * (H/2 - y))
    denominator = sp.sinh(k * H)
    psi_expression = zeta_1 * numerator / denominator

    # Create and print the equation psi(y) = expression
    final_equation = sp.Eq(psi, psi_expression)
    
    print("The derived expression for the Electrical double-layer potential distribution psi(y) is:")
    
    # Use sp.pprint for a more readable, formatted output
    sp.pprint(final_equation, use_unicode=True, wrap_line=False)

if __name__ == '__main__':
    solve_edl_potential()
