import numpy as np

# 1. Discretize the Domain and Find Parameters
delta_x = 0.3 / 3.0
delta_y = 0.3 / 3.0

lambda_val = (delta_x / delta_y)**2
alpha = lambda_val + 1

print(f"Step size Δx = {delta_x}")
print(f"Step size Δy = {delta_y}")
print(f"Coefficient λ = {lambda_val}")
print(f"Coefficient α = {alpha}\n")

# 2. Formulate the System of Equations
# The general equation is: 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)^2 * f(x,y)
# With α=2, λ=1, this becomes: 4T(i,j) - T(i,j+1) - T(i,j-1) - T(i+1,j) - T(i-1,j) = -(0.1)^2 * 100xy = -xy

# Let's map T1, T2, T3, T4 to T(i,j)
# T1 = T(i=1, j=2) at (x=0.1, y=0.2)
# T2 = T(i=2, j=2) at (x=0.2, y=0.2)
# T3 = T(i=1, j=1) at (x=0.1, y=0.1)
# T4 = T(i=2, j=1) at (x=0.2, y=0.1)

# Boundary conditions
# T(x, 0.3) = 0.5 (top, j=3)
# T(x, 0) = 0 (bottom, j=0)
# T(0.3, y) = 1 (right, i=3)
# T(0, y) = 0 (left, i=0)

# Equation for T1 (i=1, j=2): 4*T1 - T(1,3) - T3 - T2 - T(0,2) = -(0.1)*(0.2)
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.02  =>  4*T1 - T2 - T3 = 0.48

# Equation for T2 (i=2, j=2): 4*T2 - T(2,3) - T4 - T(3,2) - T1 = -(0.2)*(0.2)
# 4*T2 - 0.5 - T4 - 1 - T1 = -0.04  =>  -T1 + 4*T2 - T4 = 1.46

# Equation for T3 (i=1, j=1): 4*T3 - T1 - T(1,0) - T4 - T(0,1) = -(0.1)*(0.1)
# 4*T3 - T1 - 0 - T4 - 0 = -0.01  =>  -T1 + 4*T3 - T4 = -0.01

# Equation for T4 (i=2, j=1): 4*T4 - T2 - T(2,0) - T(3,1) - T3 = -(0.2)*(0.1)
# 4*T4 - T2 - 0 - 1 - T3 = -0.02  =>  -T2 - T3 + 4*T4 = 0.98

# 3. Solve the System in Matrix Form (A * T = B)
# The system of equations for T = [T1, T2, T3, T4]:
#  4*T1 - 1*T2 - 1*T3 + 0*T4 = 0.48
# -1*T1 + 4*T2 + 0*T3 - 1*T4 = 1.46
# -1*T1 + 0*T2 + 4*T3 - 1*T4 = -0.01
#  0*T1 - 1*T2 - 1*T3 + 4*T4 = 0.98

# Define the coefficient matrix A
A = np.array([
    [4, -1, -1, 0],
    [-1, 4, 0, -1],
    [-1, 0, 4, -1],
    [0, -1, -1, 4]
])

# Define the constant vector B
B = np.array([0.48, 1.46, -0.01, 0.98])

print("The system of linear equations is A * {T} = {B}")
print("Coefficient matrix A:")
print(A)
print("\nConstant vector B:")
print(B)

# 4. Solve for T
T = np.linalg.solve(A, B)

print("\nFinal solution for the temperature vector {T} = [T1, T2, T3, T4]:")
# Use np.array2string for cleaner formatting without extra spaces
print(np.array2string(T, separator=', ', formatter={'float_kind':lambda x: "%.4f" % x}))
