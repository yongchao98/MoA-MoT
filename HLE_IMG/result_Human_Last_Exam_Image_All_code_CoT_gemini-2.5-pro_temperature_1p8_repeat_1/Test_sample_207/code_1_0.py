import sympy

def solve_rpr_inertia():
    """
    This function provides the symbolic expression for the inertia matrix element M_32
    of the RPR robot described in the problem.
    """
    # Define the symbolic variables used in the expression.
    # m_3: mass of link 3
    # d_c3: distance to the center of mass of link 3 along its x-axis
    # q_3: joint variable for the third (revolute) joint
    m_3 = sympy.Symbol('m_3')
    d_c3 = sympy.Symbol('d_c3')
    q_3 = sympy.Symbol('q_3')
    sin = sympy.sin

    # The expression for M_32 is derived from the Lagrangian dynamics of the robot.
    # It arises from the coupling between the motion of the prismatic joint 2 and
    # the revolute joint 3 through the inertia of link 3.
    # M_32 represents the torque at joint 3 caused by a unit acceleration at joint 2.
    # Symmetrically, M_23 represents the force at joint 2 caused by a unit angular
    # acceleration at joint 3.
    M_32 = -m_3 * d_c3 * sin(q_3)

    # Print the final expression for M_32.
    # The equation is constructed by concatenating strings and the symbolic expression.
    # Each part of the equation is printed explicitly as requested.
    print("The expression for the entry M_32 of the robot inertia matrix M(q) is:")
    print(f"M_32 = {M_32.args[0]} * {M_32.args[1]} * {M_32.args[2]} * {M_32.args[3]}")
    # A more readable version:
    # print(f"M_32 = {sympy.pretty(M_32)}")
    
solve_rpr_inertia()

# Final answer to be extracted
print("<<<M_32 = -m_3*d_c3*sin(q_3)>>>")
