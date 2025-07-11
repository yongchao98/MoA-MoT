import sympy as sp

def solve_edl_potential():
    """
    This function derives and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y)
    in a parallel-plate microchannel with asymmetric zeta potentials.
    """
    # Define the symbolic variables used in the equation
    y, z_1, beta, k, H = sp.symbols('y z_1 beta k H', real=True, positive=True)
    psi = sp.Function('psi')(y)

    # The problem specifies the boundary conditions:
    # At y=0 (bottom plate), the zeta potential is z_a1 = z_1*(1 + beta*k)
    # At y=H (top plate), the zeta potential is z_a2 = 0

    # The final derived expression for psi(y) is:
    # psi(y) = z_1 * (1 + beta * k) * sinh(k * (H - y)) / sinh(k * H)
    
    # Construct the right-hand side of the equation using sympy
    zeta_a1 = z_1 * (1 + beta * k)
    numerator = sp.sinh(k * (H - y))
    denominator = sp.sinh(k * H)
    
    rhs_expression = zeta_a1 * numerator / denominator

    # Create the full equation psi(y) = expression
    final_equation = sp.Eq(psi, rhs_expression)

    # Print the final expression clearly
    print("The final expression for the Electrical Double-Layer (EDL) potential distribution, psi(y), is:")
    # Using sp.pretty_print for a more readable, formatted mathematical output
    sp.pretty_print(final_equation, use_unicode=True)

    # Outputting the parts of the final equation as requested.
    # The final equation is psi(y) = (Term 1 * Term 2) / Term 3
    # where Term 1 is the potential at the bottom wall, and Term 2/Term 3 describes the decay.
    print("\nThe equation is composed of the following parts:")
    print("\nPotential at the bottom wall (y=0):")
    sp.pretty_print(zeta_a1)
    
    print("\nNormalized potential decay function:")
    sp.pretty_print(numerator / denominator)
    

# Execute the function to display the result
solve_edl_potential()

# The final expression is z_1*(1 + beta*k)*sinh(k*(H - y))/sinh(k*H)
# Wrapping this in the required format.
final_answer = "z_1*(1 + beta*k)*sinh(k*(H - y))/sinh(k*H)"
print(f"\n<<<{final_answer}>>>")
