import math

def print_inertia_matrix_element_expression():
    """
    Prints the symbolic expression for the M_32 element of the RPR robot's inertia matrix.
    """
    
    # The derived expression for M_32.
    # It is a function of the mass of link 3 (m_3), the distance to its center of mass (d_c3),
    # and the joint angle q_3.
    expression = "-m_3 * d_c3 * sin(q_3)"
    
    print("The expression for the entry M_32 of the robot inertia matrix M(q) is:")
    print(f"M_32 = {expression}")
    
    print("\nThis expression depends on the following physical and configuration parameters:")
    print("  - m_3: The mass of link 3.")
    print("  - d_c3: The distance from the axis of joint 3 to the center of mass of link 3.")
    print("  - q_3: The angle of the third (revolute) joint.")
    
    print("\nTo represent the final equation by its components:")
    print("Each number/variable in the final equation `M_32 = -1 * m_3 * d_c3 * sin(q_3)` is:")
    print("  - Numerical coefficient: -1")
    print("  - Mass term: m_3")
    print("  - Distance term: d_c3")
    print("  - Configuration-dependent term: sin(q_3)")

if __name__ == '__main__':
    print_inertia_matrix_element_expression()
