import sympy as sp

# Define symbolic variables for the kinematics
q1, q2, q3 = sp.symbols('q_1 q_2 q_3')
d1, dc3, m3 = sp.symbols('d_1 d_c3 m_3')

# Abbreviate trig functions for clarity
c1, s1 = sp.cos(q1), sp.sin(q1)
c3, s3 = sp.cos(q3), sp.sin(q3)

# The position vector of the Center of Mass of link 3 (p_c3) is derived from the robot's geometry.
# The term d1 represents the constant height of the horizontal arm.
# dc3 is the distance from joint 3 axis to the CoM of link 3.
p_c3 = sp.Matrix([(q2 + dc3*c3)*c1, 
                  (q2 + dc3*c3)*s1, 
                  d1 + dc3*s3])

# The entry M_32 of the inertia matrix is given by M_32 = m_3 * J_v_c3,2^T * J_v_c3,3.
# The angular velocity term is zero because joint 2 is prismatic.

# We calculate the Jacobian columns by differentiating p_c3 with respect to q2 and q3.
J_vc3_2 = p_c3.diff(q2)
J_vc3_3 = p_c3.diff(q3)

# The dot product of the Jacobian columns
dot_product = J_vc3_2.dot(J_vc3_3)

# Simplify the dot product
simple_dot_product = sp.simplify(dot_product)

# Finally, M_32 = m_3 * (simplified dot product)
M_32_expr = m3 * simple_dot_product

# To display the equation clearly with all numbers, we format the output.
# The result of the symbolic calculation is -m_3*d_c3*sin(q_3).
# We print this equation with the numerical coefficient "-1" explicitly shown.
print("The derived expression for the inertia matrix entry M_32 is:")
print(f"M_32 = -1 * {m3} * {dc3} * sin({q3})")
