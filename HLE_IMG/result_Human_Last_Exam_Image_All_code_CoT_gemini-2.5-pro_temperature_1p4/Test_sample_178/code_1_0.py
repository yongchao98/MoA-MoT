import numpy as np

# 1. Define problem parameters
delta_x = 0.1
delta_y = 0.1

# 2. Calculate and print coefficients lambda and alpha
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print(f"Calculated lambda (λ): {lambda_val}")
print(f"Calculated alpha (α): {alpha_val}\n")

# 3. Define the system of linear equations AT = b
# The equations are derived from the finite difference method applied at each node.
# T_vector = [T1, T2, T3, T4]

# Coefficient matrix A
A = np.array([
    [4.0, -1.0, -1.0,  0.0],
    [-1.0,  4.0,  0.0, -1.0],
    [-1.0,  0.0,  4.0, -1.0],
    [0.0, -1.0, -1.0,  4.0]
])

# Constant vector b
b = np.array([
    0.48,
    1.46,
    -0.01,
    0.98
])

print("The system of equations is AT = b, where:")
print("\nMatrix A:")
print(A)
print("\nVector b:")
# print b as a column vector for clarity
print(b.reshape(-1, 1))

# 4. Solve the system for T
try:
    T = np.linalg.solve(A, b)

    # 5. Print the solution vector {T}
    print("\nSolving for T = [T1, T2, T3, T4]ᵀ...")
    print("\nThe solution for the temperature vector {T} is:")
    print(f"T1 = {T[0]}")
    print(f"T2 = {T[1]}")
    print(f"T3 = {T[2]}")
    print(f"T4 = {T[3]}")
    
    print("\nFinal vector value for {T} in [T1, T2, T3, T4] format:")
    print(list(T))

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")

# Final answer in the required format
print(f'<<<{[T[0], T[1], T[2], T[3]]}>>>')
