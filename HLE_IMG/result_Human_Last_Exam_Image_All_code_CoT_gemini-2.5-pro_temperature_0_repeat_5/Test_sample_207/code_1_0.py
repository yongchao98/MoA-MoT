import sympy

def solve_robot_inertia():
    """
    This function calculates and prints the symbolic expression for the M_32 element
    of the inertia matrix for the given RPR robot.
    """
    # Define the symbolic variables based on the problem description.
    # m_3: mass of link 3
    # l_3: length of link 3 (distance from joint 3 axis to joint 4 axis)
    # d_c3: distance of the Center of Mass (CoM) of link 3 from the axis of joint 4
    # q_3: joint variable for joint 3 (a revolute joint)
    m3 = sympy.Symbol('m_3')
    l3 = sympy.Symbol('l_3')
    dc3 = sympy.Symbol('d_c3')
    q3 = sympy.Symbol('q_3')

    # The distance of the CoM of link 3 from the axis of joint 3 is denoted as r_c3.
    # From the problem definition, the total length of link 3 is l_3.
    # The CoM is located at a distance d_c3 from the end of the link (joint 4 axis).
    # Therefore, its distance from the start of the link (joint 3 axis) is r_c3 = l_3 - d_c3.
    r_c3 = l3 - dc3

    # The expression for the inertia matrix element M_32 is derived from the formula:
    # M_32 = m_3 * (J_vc3_3)^T * J_vc3_2
    # After derivation, this simplifies to:
    # M_32 = -m_3 * r_c3 * cos(q_3)
    M32 = -m3 * r_c3 * sympy.cos(q3)

    # Print the final expression in a clear, readable format.
    # The prompt asks to output each number in the final equation.
    # We will explicitly show the -1 coefficient.
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    print(f"M_32 = -1 * {m3} * ({l3} - {dc3}) * cos({q3})")
    
    # For a more standard mathematical representation:
    # print(f"M_32 = {M32}")


if __name__ == "__main__":
    solve_robot_inertia()
