import numpy as np

# Step 1 & 2: Determine parameters and calculate coefficients
# The plate is a 0.3 x 0.3 square with 3 intervals in each direction.
# This means there are 4 grid points along each axis.
num_intervals = 3
domain_size = 0.3
dx = domain_size / num_intervals
dy = domain_size / num_intervals

lambda_val = (dx / dy)**2
alpha_val = lambda_val + 1

print(f"Step size dx = {dx}")
print(f"Step size dy = {dy}")
print(f"Coefficient lambda = {lambda_val}")
print(f"Coefficient alpha = {alpha_val}\n")

# Step 3 & 4: Formulate the linear equations
# The general equation is:
# 2*alpha*T(i,j) - lambda*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -dx^2 * f(x(i), y(j))
# With lambda=1 and alpha=2, this simplifies to:
# 4*T(i,j) - T(i,j+1) - T(i,j-1) - T(i+1,j) - T(i-1,j) = -dx^2 * 100 * x(i) * y(j)

# Coordinates of the unknown points:
# T1 at (x1, y2) = (0.1, 0.2)
# T2 at (x2, y2) = (0.2, 0.2)
# T3 at (x1, y1) = (0.1, 0.1)
# T4 at (x2, y1) = (0.2, 0.1)

# Let T_vec = [T1, T2, T3, T4]
# Matrix A represents the coefficients of the unknowns.
A = np.array([
    [4, -1, -1,  0],  # Equation for T1
    [-1, 4,  0, -1],  # Equation for T2
    [-1,  0, 4, -1],  # Equation for T3
    [0, -1, -1,  4]   # Equation for T4
])

# Vector b represents the constant terms from boundary conditions and the source function f(x,y).
# f(x,y) = 100xy
# RHS = -dx^2 * f(x,y) + (sum of boundary temperature terms)

# Boundary values:
T_left = 0       # T(0,y) = 0
T_right = 1      # T(0.3,y) = 1
T_bottom = 0     # T(x,0) = 0
T_top = 0.5      # T(x,0.3) = 0.5

# b vector calculation:
b1 = -(dx**2) * 100 * (0.1 * 0.2) + T_left + T_top      # For T1(0.1, 0.2), neighbors are T_left(0,0.2) and T_top(0.1,0.3)
b2 = -(dx**2) * 100 * (0.2 * 0.2) + T_right + T_top     # For T2(0.2, 0.2), neighbors are T_right(0.3,0.2) and T_top(0.2,0.3)
b3 = -(dx**2) * 100 * (0.1 * 0.1) + T_left + T_bottom   # For T3(0.1, 0.1), neighbors are T_left(0,0.1) and T_bottom(0.1,0)
b4 = -(dx**2) * 100 * (0.2 * 0.1) + T_right + T_bottom  # For T4(0.2, 0.1), neighbors are T_right(0.3,0.1) and T_bottom(0.2,0)

b = np.array([b1, b2, b3, b4])

print("The system of linear equations is A * {T} = {b}, where {T} = [T1, T2, T3, T4]")
print("Matrix A:")
print(A)
print("\nVector b:")
print(b)
print("\nExpanded form of the equations:")
print(f"Eq 1: {A[0,0]:.1f}*T1 + ({A[0,1]:.1f})*T2 + ({A[0,2]:.1f})*T3 + {A[0,3]:.1f}*T4 = {b[0]:.2f}")
print(f"Eq 2: ({A[1,0]:.1f})*T1 + {A[1,1]:.1f}*T2 + {A[1,2]:.1f}*T3 + ({A[1,3]:.1f})*T4 = {b[1]:.2f}")
print(f"Eq 3: ({A[2,0]:.1f})*T1 + {A[2,1]:.1f}*T2 + {A[2,2]:.1f}*T3 + ({A[2,3]:.1f})*T4 = {b[2]:.2f}")
print(f"Eq 4: {A[3,0]:.1f}*T1 + ({A[3,1]:.1f})*T2 + ({A[3,2]:.1f})*T3 + {A[3,3]:.1f}*T4 = {b[3]:.2f}")


# Step 5: Solve the System of Equations
T = np.linalg.solve(A, b)

# Step 6: Output the Results
print("\nSolution vector {T} = [T1, T2, T3, T4]:")
print(T)

# Format the final answer as requested
final_answer = f"[{T[0]:.5f}, {T[1]:.5f}, {T[2]:.5f}, {T[3]:.5f}]"
print(f"\nFinal Answer in the requested format:\n{final_answer}")

# Output the answer for automated checking
print(f"<<<[{T[0]},{T[1]},{T[2]},{T[3]}]>>>")