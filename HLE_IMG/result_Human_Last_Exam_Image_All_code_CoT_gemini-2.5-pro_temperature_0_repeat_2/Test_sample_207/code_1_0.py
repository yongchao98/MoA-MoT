import sympy
from sympy import sin, cos, Matrix, symbols, simplify

def solve_robotics_problem():
    """
    This function calculates the M_32 element of the inertia matrix for the RPR robot.
    """
    # Define symbolic variables
    q1, q2, q3 = symbols('q1 q2 q3')
    d1, d2, d_c3, m3 = symbols('d1 d2 d_c3 m3')

    # Helper variables for trigonometric functions
    s1, c1 = sin(q1), cos(q1)
    s3, c3 = sin(q3), cos(q3)

    # Step 1: Define the transformation matrices based on the figure's frames.
    # Let's assume the frame definition from the text which is a common interpretation.
    # Transformation from frame 1 to 0
    R1_0 = Matrix([[0, c1, -s1], [0, s1, c1], [1, 0, 0]])
    p1_0 = Matrix([0, 0, d1])

    # Transformation from frame 2 to 1
    R2_1 = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    p2_1 = Matrix([0, 0, q2])

    # Transformation from frame 3 to 2
    R3_2 = Matrix([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
    p3_2 = Matrix([d2, 0, 0])

    # Step 2: Compute the forward kinematics for the CoM of link 3.
    # Position of CoM of link 3 in frame 3 (along x3 axis)
    pc3_3 = Matrix([d_c3, 0, 0])

    # Find homogeneous transformation matrices
    A1_0 = R1_0.row_join(p1_0).col_join(Matrix([[0, 0, 0, 1]]))
    A2_1 = R2_1.row_join(p2_1).col_join(Matrix([[0, 0, 0, 1]]))
    A3_2 = R3_2.row_join(p3_2).col_join(Matrix([[0, 0, 0, 1]]))
    
    # Total transformation to frame 3
    A3_0 = A1_0 * A2_1 * A3_2
    
    # Position of CoM of link 3 in base frame
    pc3_0_h = A3_0 * pc3_3.col_join(Matrix([1]))
    pc3_0 = pc3_0_h[:3,0]

    # Step 3: Calculate the 2nd and 3rd columns of the linear velocity Jacobian for p_c3.
    Jv_c3_col2 = pc3_0.diff(q2)
    Jv_c3_col3 = pc3_0.diff(q3)

    # Step 4: Calculate the translational part of M_32.
    # M_32 = m3 * Jv_c3_col2^T * Jv_c3_col3
    dot_product = Jv_c3_col2.dot(Jv_c3_col3)
    M32_expr = m3 * dot_product
    M32_simplified = simplify(M32_expr)

    # Step 5: Print the results in a clear, step-by-step manner.
    print("The entry M_32 of the inertia matrix is given by:")
    print("M_32 = m3 * J_v_c3(:,2)^T * J_v_c3(:,3)")
    print("\nFirst, we find the position of the center of mass of link 3, p_c3^0:")
    print(f"p_c3^0 = {pc3_0}")
    print("\nNext, we find the Jacobian columns by differentiating p_c3^0 w.r.t. q2 and q3:")
    print(f"J_v_c3(:,2) = dp_c3^0/dq2 = {Jv_c3_col2}")
    print(f"J_v_c3(:,3) = dp_c3^0/dq3 = {Jv_c3_col3}")
    print("\nNow, we compute the dot product J_v_c3(:,2)^T * J_v_c3(:,3):")
    
    term1 = f"({Jv_c3_col2[0]})*({Jv_c3_col3[0]})"
    term2 = f"({Jv_c3_col2[1]})*({Jv_c3_col3[1]})"
    term3 = f"({Jv_c3_col2[2]})*({Jv_c3_col3[2]})"
    
    print(f"Dot Product = {term1} + {term2} + {term3}")
    
    simplified_term1 = simplify(Jv_c3_col2[0] * Jv_c3_col3[0])
    simplified_term2 = simplify(Jv_c3_col2[1] * Jv_c3_col3[1])
    
    print(f"            = ({simplified_term1}) + ({simplified_term2}) + (0)")
    print(f"            = {simplify(dot_product)}")
    
    print("\nFinally, we find the expression for M_32:")
    print(f"M_32 = m3 * ({simplify(dot_product)})")
    print(f"M_32 = {M32_simplified}")

solve_robotics_problem()
<<<0>>>