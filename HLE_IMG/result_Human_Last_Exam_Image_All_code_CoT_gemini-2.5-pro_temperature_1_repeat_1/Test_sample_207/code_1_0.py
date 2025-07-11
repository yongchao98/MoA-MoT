import sympy as sp

def solve_robotics_inertia():
    """
    This function defines the symbols for the robot parameters and joint variables,
    constructs the expression for the inertia matrix entry M_32, and prints it.
    """
    # Define symbolic variables for the parameters and joint variables
    # m_3: mass of link 3
    # a_2: DH parameter (length of link 3)
    # d_c3: distance to the center of mass of link 3 along its x-axis
    # q_3: joint angle of the third joint
    m_3, a_2, d_c3, q_3 = sp.symbols('m_3 a_2 d_c3 q_3')

    # Define the symbol for the inertia matrix entry
    M_32_symbol = sp.Symbol('M_32')

    # Construct the expression for M_32
    # Note: c_3 is represented as sp.cos(q_3)
    expression_M_32 = -m_3 * (a_2 + d_c3 * sp.cos(q_3))

    # Create the equation object for pretty printing
    equation = sp.Eq(M_32_symbol, expression_M_32)

    # Print the final expression
    print("The expression for the entry M_32 of the robot inertia matrix M(q) is:")
    sp.pprint(equation, use_unicode=True)

if __name__ == '__main__':
    solve_robotics_inertia()
