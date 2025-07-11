import sympy

def solve_potential_distribution():
    """
    This function defines and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y)
    based on the problem's conditions.
    """
    # Define the symbols used in the equation
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H')
    psi = sympy.Function('psi')

    # This is the derived expression for psi(y)
    # psi(y) = z_1 * (1 + beta * k) * sinh(k*(H-y)) / sinh(k*H)
    potential_expression = z_1 * (1 + beta * k) * (sympy.sinh(k * (H - y)) / sympy.sinh(k * H))

    # Create the equation psi(y) = expression
    final_equation = sympy.Eq(psi(y), potential_expression)

    # Print the final equation in a readable format
    print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
    
    # The following line prints the equation with all its symbolic components ("numbers")
    print(final_equation)

if __name__ == "__main__":
    solve_potential_distribution()