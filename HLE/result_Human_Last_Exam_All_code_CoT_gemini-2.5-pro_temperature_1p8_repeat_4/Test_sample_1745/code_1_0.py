import sympy

def solve_edl_potential():
    """
    This function defines the symbolic variables and prints the final expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y).
    """
    # Define the symbolic variables
    y, H, k, beta, z1 = sympy.symbols('y H k beta z1')
    psi = sympy.Function('psi')(y)

    # Slip-dependent zeta potential at the bottom surface (j=1)
    z_a1 = z1 * (1 + beta * k)
    
    # Final expression for the EDL potential distribution
    # psi(y) = z_a1 * sinh(k*(H-y)) / sinh(k*H)
    expression = z_a1 * sympy.sinh(k * (H - y)) / sympy.sinh(k * H)

    # Create the equation to be printed
    final_equation = sympy.Eq(psi, expression)

    # Print the final result
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    print(final_equation)
    
    # Also print the formula in the required format
    # Using sympy.sstr for a more direct string representation
    final_string_expression = sympy.sstr(final_equation, full_prec=False)
    # The output format is specified as "<<<answer content>>>"
    print(f"\n<<<{final_string_expression}>>>")

if __name__ == "__main__":
    solve_edl_potential()