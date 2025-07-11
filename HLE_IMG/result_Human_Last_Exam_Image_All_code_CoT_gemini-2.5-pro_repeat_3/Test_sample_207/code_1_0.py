import sympy as sp

def solve_robotics_problem():
    """
    This function calculates the M_32 entry of the inertia matrix for the RPR robot.
    It defines symbolic variables and prints the final expression for M_32.
    """
    # Define symbolic variables for joint angles, link parameters, and mass
    q1, q3 = sp.symbols('q1 q3')
    m3, d_c3 = sp.symbols('m3 d_c3')

    # Abbreviate trigonometric functions
    c1 = sp.cos(q1)
    s1 = sp.sin(q1)
    c3 = sp.cos(q3)
    s3 = sp.sin(q3)

    # Jacobian columns derived from the kinematics
    # J_vc3_2 represents the velocity of link 3's CoM due to joint 2's motion.
    J_vc3_2 = sp.Matrix([c1, s1, 0])
    
    # J_vc3_3 represents the velocity of link 3's CoM due to joint 3's motion.
    J_vc3_3 = sp.Matrix([d_c3 * c1 * c3, d_c3 * s1 * c3, -d_c3 * s3])

    # The inertia matrix M is symmetric, so M_32 = M_23.
    # We calculate M_23 = m3 * J_vc3_2^T * J_vc3_3
    # The angular velocity term is zero because joint 2 is prismatic.
    M_23 = m3 * (J_vc3_2.T * J_vc3_3)[0]
    
    M_32 = sp.simplify(M_23)

    # Print the final expression
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    # The sp.pretty_print function is nice but we will format it ourselves
    # to show the multiplication explicitly as requested.
    
    # Extract terms for pretty printing
    terms = M_32.as_ordered_terms()
    
    final_expr = ""
    for term in terms:
        # Separate coefficient and variables
        coeff, factors = term.as_coeff_Mul()
        
        # Build the string for factors
        factor_str = " * ".join(map(str, factors.as_ordered_factors()))
        
        if coeff == 1:
            final_expr += factor_str
        elif coeff == -1:
            final_expr += "- " + factor_str
        else:
            final_expr += f"{coeff} * {factor_str}"

    print(f"M_32 = {final_expr}")

solve_robotics_problem()