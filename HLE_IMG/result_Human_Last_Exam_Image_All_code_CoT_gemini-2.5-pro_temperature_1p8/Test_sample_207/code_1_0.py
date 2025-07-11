import sympy as sp

# This script derives the expression for the inertia matrix element M_32.
# We use sympy for symbolic mathematics to make the derivation clear.

# Define symbolic variables for the robot's parameters and joint variables.
# m_3: mass of link 3
# d_c3: distance of the Center of Mass (CoM) of link 3 along its x-axis
# q_1, q_3: joint angles for joint 1 and 3
# q_2: joint displacement for joint 2
m_3, d_c3, q_1, q_2, q_3 = sp.symbols('m_3 d_c3 q_1 q_2 q_3')

# To make the output cleaner, let's use c1, s1 for cos(q1), sin(q1) etc.
c_1, s_1 = sp.cos(q_1), sp.sin(q_1)
c_3, s_3 = sp.cos(q_3), sp.sin(q_3)

# The general formula for an inertia matrix element M_kj is:
# M_kj = Sum over i=max(k,j) to n of [ m_i * J_vi(k)^T * J_vi(j) + J_wi(k)^T * R_i * I_ci * R_i^T * J_wi(j) ]
# For M_32, this simplifies to the term for i=3:
# M_32 = m_3 * J_v3(3)^T * J_v3(2) + J_w3(3)^T * R_3 * I_c3 * R_3^T * J_w3(2)

# Step 1: Analyze the rotational component.
# The angular velocity of link 3 (omega_3) is the sum of contributions from revolute joints 1 and 3.
# omega_3 = (angular velocity from joint 1) + (angular velocity from joint 3)
# The prismatic joint 2 does not contribute to the angular velocity.
# Therefore, the angular velocity jacobian column for joint 2, J_w3(2), is the zero vector.
# This means the entire rotational component of M_32 is zero.
# M_32 = m_3 * J_v3(3)^T * J_v3(2)

# Step 2: Determine the linear Jacobian columns J_v3(2) and J_v3(3).
# We establish a kinematic model based on the figure to find the velocity of the CoM of link 3.
# Position of CoM of link 3, p_c3, is:
# p_c3 = [ (q_2 + d_c3*cos(q_3))*cos(q_1),
#          (q_2 + d_c3*cos(q_3))*sin(q_1),
#           d_1 - d_c3*sin(q_3) ]
# where d_1 is the constant height of the vertical link.

# The linear velocity v_c3 is the time derivative of p_c3.
# v_c3 = J_v3 * [q1_dot, q2_dot, q3_dot]^T
# The columns of the Jacobian J_v3 are the partial derivatives of p_c3 w.r.t q1, q2, q3.
# A more formal derivation yields the velocity coefficients for q2_dot and q3_dot:
J_v3_2_str = "[cos(q_1), sin(q_1), 0]^T"
J_v3_3_str = "[-d_c3*sin(q_3)*cos(q_1), -d_c3*sin(q_3)*sin(q_1), -d_c3*cos(q_3)]^T"

# Represent these vectors using sympy
J_v3_2 = sp.Matrix([c_1, s_1, 0])
J_v3_3 = sp.Matrix([-d_c3 * s_3 * c_1, -d_c3 * s_3 * s_1, -d_c3 * c_3])

# Step 3: Compute the dot product to find M_32.
# M_32 = m_3 * J_v3(3)^T * J_v3(2)
M_32 = m_3 * (J_v3_3.T * J_v3_2)

# The result is a 1x1 matrix, so we extract the element.
M_32_simplified = sp.simplify(M_32[0])

# Step 4: Print the derivation and the final result.
print("Derivation of the inertia matrix entry M_32:")
print("-" * 50)
print("The rotational part of M_32 is zero because joint 2 is prismatic.")
print("The formula simplifies to: M_32 = m_3 * (J_v3(3))^T * J_v3(2)\n")
print(f"The Jacobian for CoM of link 3 w.r.t. q_2 is:")
print(f"J_v3(2) = {J_v3_2_str}\n")
print(f"The Jacobian for CoM of link 3 w.r.t. q_3 is:")
print(f"J_v3(3) = {J_v3_3_str}\n")
print("Computing the dot product m_3 * (J_v3(3))^T * J_v3(2):")
print(f"M_32 = m_3 * ( (-d_c3*sin(q_3)*cos(q_1)) * cos(q_1) + (-d_c3*sin(q_3)*sin(q_1)) * sin(q_1) + (-d_c3*cos(q_3)) * 0 )")
print(f"M_32 = m_3 * ( -d_c3*sin(q_3) * (cos(q_1)^2 + sin(q_1)^2) )")
print(f"Using the identity cos(a)^2 + sin(a)^2 = 1:")
print(f"M_32 = m_3 * ( -d_c3*sin(q_3) )")
print("-" * 50)
print("\nFinal Result:")
# Use the symbolic result from sympy for the final print statement.
print(f"M_32 = {M_32_simplified}")
