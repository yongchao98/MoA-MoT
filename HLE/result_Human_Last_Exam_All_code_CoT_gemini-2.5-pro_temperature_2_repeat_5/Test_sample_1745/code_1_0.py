import sympy

def display_potential_distribution():
    """
    This function defines and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y).
    """
    # Define the symbolic variables
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H')
    psi = sympy.Function('psi')(y)

    # The problem specifies the slip-dependant zeta potential at the bottom surface (j=1)
    # z_aj = z_j(1 + beta * k)
    # The potential at the top surface is zero (z_2 = 0).
    z_a1 = z_1 * (1 + beta * k)
    
    # The derived expression for the potential distribution psi(y) is:
    # psi(y) = z_a1 * sinh(k*(H-y)) / sinh(k*H)
    expression = z_a1 * sympy.sinh(k * (H - y)) / sympy.sinh(k * H)
    
    # Create and print the final equation
    final_equation = sympy.Eq(psi, expression)
    
    # To fulfill the requirement of "output each number in the final equation",
    # we explicitly represent the numbers '1' in the equation.
    # We substitute z_a1 back into its full form for clarity.
    
    final_expression = z_1 * (sympy.Integer(1) + beta * k) * sympy.sinh(k * (H - y)) / sympy.sinh(k * H)
    final_equation_explicit = sympy.Eq(psi, final_expression)
    
    print("The expression for the Electrical double-layer potential distribution is:")
    print(sympy.pretty(final_equation_explicit, use_unicode=False))


if __name__ == '__main__':
    display_potential_distribution()
    # The expression derived is psi(y) = z_1*(1 + beta*k) * sinh(k*(H - y)) / sinh(k*H)
    # To format for the final answer block: psi(y) = z_1(1+beta*k) * sinh(k*(H-y))/sinh(kH)
    final_answer_string = "z_1(1+beta*k) * sinh(k*(H-y))/sinh(kH)"
    print(f"\n<<<{final_answer_string}>>>")