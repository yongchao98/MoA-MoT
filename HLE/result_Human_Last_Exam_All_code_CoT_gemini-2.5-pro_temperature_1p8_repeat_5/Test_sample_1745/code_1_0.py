import sympy

def solve_potential_distribution():
    """
    This function defines the symbolic variables and prints the final expression for the
    electrical double-layer potential distribution psi(y).
    """
    # Define the symbolic variables
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H')
    psi = sympy.Function('psi')(y)

    # From the problem description, the boundary conditions are:
    # 1. psi(y = H/2) = 0
    # 2. psi(y = -H/2) = z_1 * (1 + beta*k)
    # The solution to the linearized Poisson-Boltzmann equation d^2(psi)/dy^2 = k^2 * psi
    # with these boundary conditions is derived as follows:
    
    # The numerator of the expression
    numerator = z_1 * (1 + beta * k) * sympy.sinh(k * (H/2 - y))
    
    # The denominator of the expression
    denominator = sympy.sinh(k * H)
    
    # The final expression for psi(y)
    psi_expression = numerator / denominator
    
    # Create and print the equation object
    # The output shows all symbolic variables in the final derived equation
    final_equation = sympy.Eq(psi, psi_expression)
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    print(final_equation)

solve_potential_distribution()