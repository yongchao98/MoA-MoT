def solve_rpr_inertia():
    """
    This function provides the derived expression for the M_32 element
    of the inertia matrix for the given RPR spatial robot.
    """
    # The derived symbolic expression is M_32 = -m_3 * d_c3 * sin(q_3).
    # We will print this expression in a clear format.
    # The variable m_3 is the mass of link 3.
    # The variable d_c3 is the distance of the Center of Mass of link 3 along the x3 axis.
    # The variable q_3 is the joint angle for joint 3.
    
    print("The expression for the entry M_32 of the robot inertia matrix M(q) is:")
    
    # We represent the expression as a string.
    # To satisfy the instruction "output each number in the final equation",
    # we explicitly write the negative sign as "-1 *".
    expression = "-1 * m_3 * d_c3 * sin(q_3)"
    
    print("M_32 =", expression)

solve_rpr_inertia()