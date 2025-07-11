import numpy as np

# Step 1 & 2: Define parameters and calculate coefficients
delta_x = 0.1
delta_y = 0.1
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print("Step 1: Calculated Coefficients")
print(f"λ = ({delta_x}/{delta_y})² = {lambda_val}")
print(f"α = λ + 1 = {lambda_val} + 1 = {alpha_val}\n")

# Step 3 & 4: Formulate the system of linear equations A * {T} = {b}
# The general equation is: 4*T_center - T_north - T_south - T_east - T_west = -x*y
# The vector of unknowns is {T} = [T1, T2, T3, T4]

# Boundary values: T(left)=0, T(right)=1, T(bottom)=0, T(top)=0.5

# Equation for T1 at (x=0.1, y=0.2):
# 4*T1 - T(top) - T3 - T2 - T(left) = -(0.1)*(0.2)
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.02
# 4*T1 - T2 - T3 = 0.48
eq1_coeffs = [4., -1., -1., 0.]
eq1_b = 0.48

# Equation for T2 at (x=0.2, y=0.2):
# 4*T2 - T(top) - T4 - T(right) - T1 = -(0.2)*(0.2)
# 4*T2 - 0.5 - T4 - 1.0 - T1 = -0.04
# -T1 + 4*T2 - T4 = 1.46
eq2_coeffs = [-1., 4., 0., -1.]
eq2_b = 1.46

# Equation for T3 at (x=0.1, y=0.1):
# 4*T3 - T1 - T(bottom) - T4 - T(left) = -(0.1)*(0.1)
# 4*T3 - T1 - 0 - T4 - 0 = -0.01
# -T1 + 4*T3 - T4 = -0.01
eq3_coeffs = [-1., 0., 4., -1.]
eq3_b = -0.01

# Equation for T4 at (x=0.2, y=0.1):
# 4*T4 - T2 - T(bottom) - T(right) - T3 = -(0.2)*(0.1)
# 4*T4 - T2 - 0 - 1.0 - T3 = -0.02
# -T2 - T3 + 4*T4 = 0.98
eq4_coeffs = [0., -1., -1., 4.]
eq4_b = 0.98

# Assemble the matrix A and vector b
A = np.array([eq1_coeffs, eq2_coeffs, eq3_coeffs, eq4_coeffs])
b = np.array([eq1_b, eq2_b, eq3_b, eq4_b])

print("Step 2: The system of linear equations in the form A*{T}={b} is:")
print(f"({A[0,0]})*T1 + ({A[0,1]})*T2 + ({A[0,2]})*T3 + ({A[0,3]})*T4 = {b[0]}")
print(f"({A[1,0]})*T1 + ({A[1,1]})*T2 + ({A[1,2]})*T3 + ({A[1,3]})*T4 = {b[1]}")
print(f"({A[2,0]})*T1 + ({A[2,1]})*T2 + ({A[2,2]})*T3 + ({A[2,3]})*T4 = {b[2]}")
print(f"({A[3,0]})*T1 + ({A[3,1]})*T2 + ({A[3,2]})*T3 + ({A[3,3]})*T4 = {b[3]}\n")

# Step 5: Solve the system A * T = b
T = np.linalg.solve(A, b)

# Step 6: Present the solution
print("Step 3: The solution vector {T} is:")
print(f"[{T[0]}, {T[1]}, {T[2]}, {T[3]}]")