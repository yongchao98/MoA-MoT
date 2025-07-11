import numpy as np

# Step 1: Define problem parameters and calculate coefficients
# The domain is a 0.3 x 0.3 square with 3 intervals in each direction.
delta_x = 0.3 / 3.0
delta_y = 0.3 / 3.0

# Calculate lambda and alpha
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print("Step 1: Calculate coefficients")
print(f"Δx = {delta_x}")
print(f"Δy = {delta_y}")
print(f"λ = (Δx/Δy)² = {lambda_val}")
print(f"α = λ + 1 = {alpha_val}\n")

# Step 2: Set up the system of linear equations A * {T} = {B}
# The unknowns are {T} = [T1, T2, T3, T4]^T.
# The finite difference equation is:
# 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)² * f(x,y)
# 4*T(i,j) - T_up - T_down - T_right - T_left = -(0.1)² * 100*x*y = -x*y

# Boundary Conditions:
# T(0,y) = 0 (left), T(0.3,y) = 1 (right)
# T(x,0) = 0 (bottom), T(x,0.3) = 0.5 (top)

# Node T1 (x=0.1, y=0.2):
# 4*T1 - T(0.1,0.3) - T(0.1,0.1) - T(0.2,0.2) - T(0,0.2) = -(0.1)*(0.2)
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.02
# 4*T1 - T2 - T3 = 0.48

# Node T2 (x=0.2, y=0.2):
# 4*T2 - T(0.2,0.3) - T(0.2,0.1) - T(0.3,0.2) - T(0.1,0.2) = -(0.2)*(0.2)
# 4*T2 - 0.5 - T4 - 1 - T1 = -0.04
# -T1 + 4*T2 - T4 = 1.46

# Node T3 (x=0.1, y=0.1):
# 4*T3 - T(0.1,0.2) - T(0.1,0) - T(0.2,0.1) - T(0,0.1) = -(0.1)*(0.1)
# 4*T3 - T1 - 0 - T4 - 0 = -0.01
# -T1 + 4*T3 - T4 = -0.01

# Node T4 (x=0.2, y=0.1):
# 4*T4 - T(0.2,0.2) - T(0.2,0) - T(0.3,0.1) - T(0.1,0.1) = -(0.2)*(0.1)
# 4*T4 - T2 - 0 - 1 - T3 = -0.02
# -T2 - T3 + 4*T4 = 0.98

# Coefficient matrix A
A = np.array([
    [ 4, -1, -1,  0],
    [-1,  4,  0, -1],
    [-1,  0,  4, -1],
    [ 0, -1, -1,  4]
])

# Constant vector B
B = np.array([0.48, 1.46, -0.01, 0.98])

print("Step 2: Formulate the system of linear equations A * {T} = {B}")
print("The unknowns are {T} = [T1, T2, T3, T4]^T")
print("The derived equations are:")
print(f"({A[0,0]:.0f})*T1 + ({A[0,1]:.0f})*T2 + ({A[0,2]:.0f})*T3 + ({A[0,3]:.0f})*T4 = {B[0]}")
print(f"({A[1,0]:.0f})*T1 + ({A[1,1]:.0f})*T2 + ({A[1,2]:.0f})*T3 + ({A[1,3]:.0f})*T4 = {B[1]}")
print(f"({A[2,0]:.0f})*T1 + ({A[2,1]:.0f})*T2 + ({A[2,2]:.0f})*T3 + ({A[2,3]:.0f})*T4 = {B[2]}")
print(f"({A[3,0]:.0f})*T1 + ({A[3,1]:.0f})*T2 + ({A[3,2]:.0f})*T3 + ({A[3,3]:.0f})*T4 = {B[3]}\n")

# Step 3: Solve the system for T
T = np.linalg.solve(A, B)

# Step 4: Print the result
print("Step 3: Solve for {T}")
print(f"The vector value for {{T}} is: [{T[0]:.4f}, {T[1]:.4f}, {T[2]:.4f}, {T[3]:.4f}]")