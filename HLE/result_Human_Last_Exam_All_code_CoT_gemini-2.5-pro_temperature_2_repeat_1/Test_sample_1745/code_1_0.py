import sympy

def solve_edl_potential():
    """
    This script symbolically derives and prints the expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y)
    in a parallel-plate microchannel under the given conditions.
    """
    # Define the mathematical symbols involved in the equation.
    # y: spatial coordinate across the channel height
    # z_1: base zeta potential on the bottom surface
    # beta: slip length
    # k: Debye-Huckel parameter
    # H: height of the microchannel
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H', real=True, positive=True)

    # Define psi as a function of y
    psi = sympy.Function('psi')(y)

    # The problem specifies a slip-dependant zeta potential at the bottom surface (y=0)
    # as z_a1 = z_1*(1 + beta*k).
    # The zeta potential at the top surface (y=H) is zero (since z_2=0).
    # The governing equation is the linearized Poisson-Boltzmann equation: d^2(psi)/dy^2 = k^2 * psi
    # The solution, after applying the boundary conditions psi(0) = z_1*(1+beta*k)
    # and psi(H) = 0, is derived as shown in the explanation.

    # Construct the final derived expression for the potential psi(y)
    numerator = sympy.sinh(k * (H - y))
    denominator = sympy.sinh(k * H)
    potential_expression = z_1 * (1 + beta * k) * (numerator / denominator)

    # Create the final equation psi(y) = expression
    final_equation = sympy.Eq(psi, potential_expression)

    # Print the final equation
    print("The derived expression for the Electrical double-layer (EDL) potential distribution psi(y) is:")
    sympy.pprint(final_equation, use_unicode=True)
    
    # Return the string representation for the final answer format
    return str(final_equation)

# Execute the function and capture the string output
final_answer_string = solve_edl_potential()
# The final answer needs to be wrapped according to instructions.
# Let's extract just the right hand side of the equation.
rhs = final_answer_string.split("=")[1].strip()
print(f"\n<<<psi(y) = {rhs}>>>")
