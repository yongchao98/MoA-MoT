import numpy as np

# 1. Determine Grid Parameters
delta_x = 0.1
delta_y = 0.1

# 2. Calculate Coefficients
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

# 3. Formulate and 4. Solve the System of Equations
# The general finite difference equation is:
# 2*alpha*T(i,j) - lambda*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(dx)^2 * f(x(i), y(j))
#
# Substituting alpha=2, lambda=1, dx=0.1, and f(x,y)=100xy simplifies to:
# 4*T(i,j) - T_up - T_down - T_right - T_left = -x(i)*y(j)
#
# Let the unknown temperatures be T_vec = [T1, T2, T3, T4]
#
# Equation for T1 at (x=0.1, y=0.2):
# 4*T1 - T(x=0.1, y=0.3) - T3 - T2 - T(x=0, y=0.2) = -0.1*0.2
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.02  =>  4*T1 - T2 - T3 = 0.48
#
# Equation for T2 at (x=0.2, y=0.2):
# 4*T2 - T(x=0.2, y=0.3) - T4 - T(x=0.3, y=0.2) - T1 = -0.2*0.2
# 4*T2 - 0.5 - T4 - 1 - T1 = -0.04  => -T1 + 4*T2 - T4 = 1.46
#
# Equation for T3 at (x=0.1, y=0.1):
# 4*T3 - T1 - T(x=0.1, y=0) - T4 - T(x=0, y=0.1) = -0.1*0.1
# 4*T3 - T1 - 0 - T4 - 0 = -0.01  => -T1 + 4*T3 - T4 = -0.01
#
# Equation for T4 at (x=0.2, y=0.1):
# 4*T4 - T2 - T(x=0.2, y=0) - T(x=0.3, y=0.1) - T3 = -0.2*0.1
# 4*T4 - T2 - 0 - 1 - T3 = -0.02  => -T2 - T3 + 4*T4 = 0.98
#
# This gives the matrix system A*T = b

# Define the coefficient matrix A
A = np.array([
    [ 4.0, -1.0, -1.0,  0.0],
    [-1.0,  4.0,  0.0, -1.0],
    [-1.0,  0.0,  4.0, -1.0],
    [ 0.0, -1.0, -1.0,  4.0]
])

# Define the constant vector b
b = np.array([0.48, 1.46, -0.01, 0.98])

# Solve the linear system A * T = b for T
T_vector = np.linalg.solve(A, b)

# Print the results
print(f"The calculated coefficients are: lambda = {lambda_val:.1f} and alpha = {alpha_val:.1f}")
print("\nThe system of linear equations derived for the nodes T1, T2, T3, and T4 is:")
print(f"  {A[0,0]:.1f}*T1 + ({A[0,1]:.1f})*T2 + ({A[0,2]:.1f})*T3 + {A[0,3]:.1f}*T4  = {b[0]:.2f}")
print(f" ({A[1,0]:.1f})*T1 + {A[1,1]:.1f}*T2 + {A[1,2]:.1f}*T3 + ({A[1,3]:.1f})*T4  = {b[1]:.2f}")
print(f" ({A[2,0]:.1f})*T1 + {A[2,1]:.1f}*T2 + {A[2,2]:.1f}*T3 + ({A[2,3]:.1f})*T4 = {b[2]:.2f}")
print(f"  {A[3,0]:.1f}*T1 + ({A[3,1]:.1f})*T2 + ({A[3,2]:.1f})*T3 + {A[3,3]:.1f}*T4  = {b[3]:.2f}")
print("\nSolving this system gives the temperature vector {T}.")
print(f"\nThe vector value for {{T}} is: [T1, T2, T3, T4] = {T_vector.tolist()}")
