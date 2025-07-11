import sympy

def solve_inertia_matrix_element():
    """
    This function derives and prints the expression for the M_32 element
    of the inertia matrix for the given RPR robot.
    """

    # Define the symbolic variables based on the problem description.
    # m_3: mass of link 3
    # d_c3: distance of the Center of Mass of link 3 from joint 3 axis, along the x3 axis
    # q_3: joint variable 3 (rotation angle)
    m_3 = sympy.Symbol('m_3')
    d_c3 = sympy.Symbol('d_c3')
    q_3 = sympy.Symbol('q_3')

    # The derived expression for M_32 is -m_3 * d_c3 * sin(q_3).
    # We construct this expression and print it.
    
    # The final expression includes the symbolic variables and a coefficient of -1.
    # The components are the coefficient -1, the variable m_3, the variable d_c3,
    # and the function sin(q_3).
    
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    print(f"M_32 = -{m_3} * {d_c3} * sin({q_3})")

solve_inertia_matrix_element()