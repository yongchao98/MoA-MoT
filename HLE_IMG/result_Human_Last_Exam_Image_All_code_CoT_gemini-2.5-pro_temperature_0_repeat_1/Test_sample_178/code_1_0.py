import numpy as np

# Step 1 & 2: Discretization and Coefficients
delta_x = 0.1
delta_y = 0.1
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print("Step 1 & 2: Discretization and Coefficients")
print(f"Δx = {delta_x}")
print(f"Δy = {delta_y}")
print(f"λ = (Δx/Δy)^2 = {lambda_val}")
print(f"α = λ + 1 = {alpha_val}\n")

# The finite difference equation is:
# 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)^2 * f(x(i), y(j))
# Substituting λ=1, α=2, Δx=0.1, and f(x,y)=100xy:
# 4*T(i,j) - T(i,j+1) - T(i,j-1) - T(i+1,j) - T(i-1,j) = -x(i)*y(j)

print("Step 3: Setting up the System of Linear Equations")
print("The general equation is: 4*T_ij - T_i,j+1 - T_i,j-1 - T_i+1,j - T_i-1,j = -x_i*y_j\n")

# Boundary Conditions
# T(0,y) = 0, T(0.3,y) = 1
# T(x,0) = 0, T(x,0.3) = 0.5

# Let T_vec = [T1, T2, T3, T4]
# T1 = T(1,2), T2 = T(2,2), T3 = T(1,1), T4 = T(2,1)

# Equation for T1 (i=1, j=2) at (x=0.1, y=0.2)
# 4*T1 - T(1,3) - T(1,1) - T(2,2) - T(0,2) = -(0.1)*(0.2)
# T(1,3) is boundary T(0.1, 0.3) = 0.5
# T(1,1) is T3
# T(2,2) is T2
# T(0,2) is boundary T(0, 0.2) = 0
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.02  =>  4*T1 - T2 - T3 = 0.48
print("Equation for T1: 4 * T1 - 1 * T2 - 1 * T3 + 0 * T4 = 0.48")

# Equation for T2 (i=2, j=2) at (x=0.2, y=0.2)
# 4*T2 - T(2,3) - T(2,1) - T(3,2) - T(1,2) = -(0.2)*(0.2)
# T(2,3) is boundary T(0.2, 0.3) = 0.5
# T(2,1) is T4
# T(3,2) is boundary T(0.3, 0.2) = 1
# T(1,2) is T1
# 4*T2 - 0.5 - T4 - 1 - T1 = -0.04  =>  -T1 + 4*T2 - T4 = 1.46
print("Equation for T2: -1 * T1 + 4 * T2 + 0 * T3 - 1 * T4 = 1.46")

# Equation for T3 (i=1, j=1) at (x=0.1, y=0.1)
# 4*T3 - T(1,2) - T(1,0) - T(2,1) - T(0,1) = -(0.1)*(0.1)
# T(1,2) is T1
# T(1,0) is boundary T(0.1, 0) = 0
# T(2,1) is T4
# T(0,1) is boundary T(0, 0.1) = 0
# 4*T3 - T1 - 0 - T4 - 0 = -0.01  =>  -T1 + 4*T3 - T4 = -0.01
print("Equation for T3: -1 * T1 + 0 * T2 + 4 * T3 - 1 * T4 = -0.01")

# Equation for T4 (i=2, j=1) at (x=0.2, y=0.1)
# 4*T4 - T(2,2) - T(2,0) - T(3,1) - T(1,1) = -(0.2)*(0.1)
# T(2,2) is T2
# T(2,0) is boundary T(0.2, 0) = 0
# T(3,1) is boundary T(0.3, 0.1) = 1
# T(1,1) is T3
# 4*T4 - T2 - 0 - 1 - T3 = -0.02  =>  -T2 - T3 + 4*T4 = 0.98
print("Equation for T4: 0 * T1 - 1 * T2 - 1 * T3 + 4 * T4 = 0.98\n")

# Step 4: Solve the System A * {T} = {B}
# Coefficient matrix A
A = np.array([
    [ 4, -1, -1,  0],
    [-1,  4,  0, -1],
    [-1,  0,  4, -1],
    [ 0, -1, -1,  4]
])

# Constant vector B
B = np.array([0.48, 1.46, -0.01, 0.98])

# Solve for T
T = np.linalg.solve(A, B)

# Step 5: Output the Result
print("Step 4 & 5: Solve and Display the Result")
print("The vector value for {T} is:")
print(f"[T1, T2, T3, T4] = {T.tolist()}")
