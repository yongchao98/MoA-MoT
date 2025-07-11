import numpy as np

# Step 1: Define grid and model parameters
delta_x = 0.1
delta_y = 0.1

lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print(f"Calculated Parameters:")
print(f"λ = {lambda_val}")
print(f"α = {alpha_val}\n")

# Step 2: Define the source function f(x, y)
def f(x, y):
    return 100 * x * y

# Step 3: Set up the system of linear equations A * T = b

# Coordinates of the unknown temperature points
coords = {
    'T1': (0.1, 0.2),
    'T2': (0.2, 0.2),
    'T3': (0.1, 0.1),
    'T4': (0.2, 0.1)
}

# Boundary conditions
# T(0, y) = 0
# T(0.3, y) = 1
# T(x, 0) = 0
# T(x, 0.3) = 0.5
T_left = 0
T_right = 1
T_bottom = 0
T_top = 0.5

# Assemble the coefficient matrix A
#      T1  T2  T3  T4
A = np.array([
    [ 4, -1, -1,  0],  # Equation for T1
    [-1,  4,  0, -1],  # Equation for T2
    [-1,  0,  4, -1],  # Equation for T3
    [ 0, -1, -1,  4]   # Equation for T4
])

# Assemble the constant vector b
# For each equation, b_i = sum(known_boundary_neighbors) - (delta_x**2) * f(x_i, y_i)
b1 = T_left + T_top - (delta_x**2) * f(coords['T1'][0], coords['T1'][1])
b2 = T_right + T_top - (delta_x**2) * f(coords['T2'][0], coords['T2'][1])
b3 = T_left + T_bottom - (delta_x**2) * f(coords['T3'][0], coords['T3'][1])
b4 = T_right + T_bottom - (delta_x**2) * f(coords['T4'][0], coords['T4'][1])
b = np.array([b1, b2, b3, b4])

print("System of linear equations (A * T = b):")
print("Matrix A:")
print(A)
print("\nVector b:")
print(b)
print("\nWhich corresponds to the following equations:")
print(f"{A[0,0]}*T1 + {A[0,1]}*T2 + {A[0,2]}*T3 + {A[0,3]}*T4 = {b[0]:.2f}")
print(f"{A[1,0]}*T1 + {A[1,1]}*T2 + {A[1,2]}*T3 + {A[1,3]}*T4 = {b[1]:.2f}")
print(f"{A[2,0]}*T1 + {A[2,1]}*T2 + {A[2,2]}*T3 + {A[2,3]}*T4 = {b[2]:.2f}")
print(f"{A[3,0]}*T1 + {A[3,1]}*T2 + {A[3,2]}*T3 + {A[3,3]}*T4 = {b[3]:.2f}\n")


# Step 4: Solve the system for T
T = np.linalg.solve(A, b)

print("Solution for the temperature vector {T}:")
print(f"[T1, T2, T3, T4] = {T.tolist()}")
