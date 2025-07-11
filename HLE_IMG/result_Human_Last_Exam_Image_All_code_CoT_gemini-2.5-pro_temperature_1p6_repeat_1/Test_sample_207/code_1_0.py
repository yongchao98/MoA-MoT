import sympy as sp

def main():
    """
    Calculates the M_32 entry of the inertia matrix for the RPR robot.
    """
    # Define symbolic variables
    q1, q3 = sp.symbols('q1 q3')
    m3, I3, dc3 = sp.symbols('m3 I3 d_c3')

    # Trigonometric functions
    c1, s1 = sp.cos(q1), sp.sin(q1)
    c3, s3 = sp.cos(q3), sp.sin(q3)

    # === Step 1: Define Joint Axes in base frame {0} ===
    # Axis for prismatic joint 2 (along x1)
    z_j2 = sp.Matrix([c1, s1, 0])
    
    # Axis for revolute joint 3 (along z2, which is -x1)
    z_j3 = -z_j2

    # === Step 2: Calculate Jacobian Columns for Link 3 ===
    
    # -- Column 2 (due to joint 2) --
    # Prismatic joint: linear velocity is along the axis, no angular velocity.
    J_vc3_col2 = z_j2
    J_omega3_col2 = sp.Matrix([0, 0, 0])

    # -- Column 3 (due to joint 3) --
    # Revolute joint: angular velocity is along the axis.
    J_omega3_col3 = z_j3
    
    # To find J_vc3_col3, we need the vector from joint 3 origin to CoM of link 3.
    # The CoM of link 3 is dc3 along the x3 axis. We need the x3 direction vector.
    # We can use the physical argument that the dot product will be zero, 
    # but for completeness, let's derive it.
    # r = p_c3 - p_j3. We showed this vector is perpendicular to z_j3.
    # Let's verify the dot product directly.
    # The dot product J_vc3_col3.dot(J_vc3_col2) is equivalent to (z_j3 x r).dot(z_j2).
    # Since z_j2 = -z_j3, this is (z_j3 x r).dot(-z_j3).
    # This is a scalar triple product where two vectors are parallel, so the result is 0.
    
    # === Step 3: Compute Dot Products for M_32 ===
    # Linear part contribution
    linear_dot_product = J_vc3_col2.dot(z_j3.cross(sp.Matrix(['rx','ry','rz']))) # The exact r vector is complex but not needed
    # As z_j3 is parallel to J_vc3_col2, the triple product is zero. Let's explicitly compute.
    # J_vc3_col2.dot(J_vc3_col3) = 0
    linear_term_val = 0 # Based on vector properties
    
    # Angular part contribution
    angular_dot_product = J_omega3_col2.dot(J_omega3_col3)
    angular_term_val = sp.simplify(angular_dot_product)

    # === Step 4: Calculate M_32 ===
    M32 = m3 * linear_term_val + I3 * angular_term_val

    # Print the result in equation form
    print("The expression for the entry M_32 of the robot inertia matrix is derived from:")
    print("M_32 = m_3 * (J_v_c3[:,2]^T * J_v_c3[:,3]) + I_3 * (J_omega_3[:,2]^T * J_omega_3[:,3])")
    print("\nBreaking down the terms:")
    print("Linear part dot product J_v_c3[:,2]^T * J_v_c3[:,3] = ", linear_term_val)
    print("Angular part dot product J_omega_3[:,2]^T * J_omega_3[:,3] = ", angular_term_val)
    print("\nSubstituting these values into the equation:")
    print(f"M_32 = m3 * ({linear_term_val}) + I3 * ({angular_term_val})")
    print(f"M_32 = {m3*linear_term_val} + {I3*angular_term_val}")
    print("\nFinal expression:")
    print(f"M_32 = {M32}")

if __name__ == '__main__':
    main()
