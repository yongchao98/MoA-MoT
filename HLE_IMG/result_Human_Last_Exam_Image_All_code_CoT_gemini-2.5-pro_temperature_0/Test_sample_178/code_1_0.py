import numpy as np

# Step 1 & 2: Define geometry, discretization, and coordinates
dx = 0.1
dy = 0.1
coords = {
    'T1': (0.1, 0.2),
    'T2': (0.2, 0.2),
    'T3': (0.1, 0.1),
    'T4': (0.2, 0.1)
}

# Step 3: Calculate and print coefficients
lambda_val = (dx / dy)**2
alpha = lambda_val + 1
print(f"Calculated coefficients:")
print(f"λ = (Δx/Δy)² = ({dx}/{dy})² = {lambda_val:.2f}")
print(f"α = λ + 1 = {lambda_val:.2f} + 1 = {alpha:.2f}\n")

# Define the source function
def f(x, y):
    return 100 * x * y

# Step 4: Formulate the system of linear equations A*T = B
# The general equation is:
# 2*alpha*T_ij - lambda*[T_up + T_down] - [T_right + T_left] = -dx^2 * f(x,y)
# With lambda=1, alpha=2, this simplifies to:
# 4*T_ij - T_up - T_down - T_right - T_left = -dx^2 * f(x,y)

# Boundary Conditions:
# T(0, y) = 0 (left)
# T(0.3, y) = 1 (right)
# T(x, 0) = 0 (bottom)
# T(x, 0.3) = 0.5 (top)

# The system is ordered T1, T2, T3, T4
# A is the coefficient matrix for [T1, T2, T3, T4]
A = np.array([
    [4., -1., -1., 0.],
    [-1., 4., 0., -1.],
    [-1., 0., 4., -1.],
    [0., -1., -1., 4.]
])

# B is the vector of constants (RHS)
# B_i = (boundary terms) - dx^2 * f(x_i, y_i)
B = np.zeros(4)

# Equation for T1 (x=0.1, y=0.2)
# 4*T1 - T_top - T3 - T2 - T_left = -dx^2*f(0.1, 0.2)
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.01 * 100 * 0.1 * 0.2
# 4*T1 - T2 - T3 = 0.5 - 0.02
B[0] = 0.5 - (dx**2 * f(coords['T1'][0], coords['T1'][1]))

# Equation for T2 (x=0.2, y=0.2)
# 4*T2 - T_top - T4 - T_right - T1 = -dx^2*f(0.2, 0.2)
# 4*T2 - 0.5 - T4 - 1.0 - T1 = -0.01 * 100 * 0.2 * 0.2
# -T1 + 4*T2 - T4 = 1.5 - 0.04
B[1] = 0.5 + 1.0 - (dx**2 * f(coords['T2'][0], coords['T2'][1]))

# Equation for T3 (x=0.1, y=0.1)
# 4*T3 - T1 - T_bottom - T4 - T_left = -dx^2*f(0.1, 0.1)
# 4*T3 - T1 - 0 - T4 - 0 = -0.01 * 100 * 0.1 * 0.1
# -T1 + 4*T3 - T4 = -0.01
B[2] = 0 + 0 - (dx**2 * f(coords['T3'][0], coords['T3'][1]))

# Equation for T4 (x=0.2, y=0.1)
# 4*T4 - T2 - T_bottom - T_right - T3 = -dx^2*f(0.2, 0.1)
# 4*T4 - T2 - 0 - 1.0 - T3 = -0.01 * 100 * 0.2 * 0.1
# -T2 - T3 + 4*T4 = 1.0 - 0.02
B[3] = 0 + 1.0 - (dx**2 * f(coords['T4'][0], coords['T4'][1]))

print("The system of linear equations is A * {T} = {B}:")
vars = ['T1', 'T2', 'T3', 'T4']
for i in range(4):
    equation_str = ""
    for j in range(4):
        equation_str += f"{A[i, j]:+6.2f} * {vars[j]} "
    equation_str += f"= {B[i]:.4f}"
    print(equation_str)

# Step 5: Solve the system
T = np.linalg.solve(A, B)

# Step 6: Output the result
print("\nSolving the system gives the vector {T}:")
print(f"The vector value for {{T}} is: [{T[0]:.5f}, {T[1]:.5f}, {T[2]:.5f}, {T[3]:.5f}]")

# Final answer in the required format
answer = f"[{T[0]:.5f}, {T[1]:.5f}, {T[2]:.5f}, {T[3]:.5f}]"
# print(f"\n<<<{answer}>>>")