import numpy as np

# 1. Define parameters and boundary conditions
dx = 0.1
dy = 0.1
lam = (dx / dy)**2
alpha = lam + 1

# Coordinates of the unknown points
coords = {
    'T1': (0.1, 0.2),
    'T2': (0.2, 0.2),
    'T3': (0.1, 0.1),
    'T4': (0.2, 0.1)
}

# Boundary Conditions
T_left = 0      # T(0, y) = 0
T_right = 1     # T(0.3, y) = 1
T_bottom = 0    # T(x, 0) = 0
T_top = 0.5     # T(x, 0.3) = 0.5

# 2. Set up the system of linear equations A*T = b
# The matrix A is based on the coefficients of T1, T2, T3, T4
A = np.array([
    [4, -1, -1,  0],  # Equation for T1
    [-1, 4,  0, -1],  # Equation for T2
    [-1, 0,  4, -1],  # Equation for T3
    [0, -1, -1,  4]   # Equation for T4
])

# The vector b is calculated from the boundary conditions and the source term f(x,y)
# f(x,y) = 100xy, so the source term contribution is -(dx**2) * 100 * x * y = -x*y
b1 = - (coords['T1'][0] * coords['T1'][1]) + T_top + T_left
b2 = - (coords['T2'][0] * coords['T2'][1]) + T_top + T_right
b3 = - (coords['T3'][0] * coords['T3'][1]) + T_bottom + T_left
b4 = - (coords['T4'][0] * coords['T4'][1]) + T_bottom + T_right

b = np.array([b1, b2, b3, b4])

print("The system of linear equations is AT = b, with T = [T1, T2, T3, T4]'")
print("The equations are:")
print(f"({A[0,0]})*T1 + ({A[0,1]})*T2 + ({A[0,2]})*T3 + ({A[0,3]})*T4 = {b[0]:.2f}")
print(f"({A[1,0]})*T1 + ({A[1,1]})*T2 + ({A[1,2]})*T3 + ({A[1,3]})*T4 = {b[1]:.2f}")
print(f"({A[2,0]})*T1 + ({A[2,1]})*T2 + ({A[2,2]})*T3 + ({A[2,3]})*T4 = {b[2]:.2f}")
print(f"({A[3,0]})*T1 + ({A[3,1]})*T2 + ({A[3,2]})*T3 + ({A[3,3]})*T4 = {b[3]:.2f}")
print("")

# 3. Solve the system for the temperature vector T
T = np.linalg.solve(A, b)

# 4. Print the final result
print("The vector value for {T} is:")
# Use np.round to avoid floating point inaccuracies in the final display
T_rounded = np.round(T, 5)
print(f"{T_rounded.tolist()}")