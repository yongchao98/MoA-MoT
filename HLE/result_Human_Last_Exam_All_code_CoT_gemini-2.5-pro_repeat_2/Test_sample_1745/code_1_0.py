import sympy

def display_potential_distribution():
    """
    This function defines and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y).
    """
    # Define the mathematical symbols
    y, H, k, beta = sympy.symbols('y H k beta', real=True, positive=True)
    zeta_1 = sympy.Symbol('zeta_1', real=True)
    psi = sympy.Function('psi')

    # The slip-dependent zeta potential at the bottom surface (j=1)
    # The number 1 is explicitly used here.
    zeta_a1 = zeta_1 * (1 + beta * k)

    # The derived expression for psi(y).
    # The number 2 is explicitly used in H/2.
    numerator = sympy.sinh(k * (H/2 - y))
    denominator = sympy.sinh(k * H)
    expression = zeta_a1 * numerator / denominator

    # Create the final equation psi(y) = expression
    final_equation = sympy.Eq(psi(y), expression)

    # Print the final equation in a readable format.
    # The sympy.pretty_print function provides a nicely formatted output.
    print("The expression for the Electrical Double-Layer (EDL) potential distribution psi(y) is:")
    sympy.pretty_print(final_equation, use_unicode=True)

if __name__ == '__main__':
    display_potential_distribution()