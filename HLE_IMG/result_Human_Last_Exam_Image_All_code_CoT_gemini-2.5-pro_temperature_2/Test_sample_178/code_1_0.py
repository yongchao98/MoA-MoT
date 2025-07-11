import numpy as np

# Step 1 & 2: Define parameters and calculate coefficients
dx = 0.1
dy = 0.1

lambda_val = (dx / dy)**2
alpha_val = lambda_val + 1

print("--- Step 1 & 2: Parameters and Coefficients ---")
print(f"Grid spacing Δx = {dx}, Δy = {dy}")
print(f"λ = (Δx/Δy)² = ({dx}/{dy})² = {lambda_val:.4f}")
print(f"α = λ + 1 = {lambda_val:.4f} + 1 = {alpha_val:.4f}\n")

# Step 3: Define source function, boundary conditions, and node coordinates
def f(x, y):
    return 100 * x * y

# Boundary Temperatures
T_left = 0      # T(0, y) = 0
T_right = 1     # T(0.3, y) = 1
T_bottom = 0    # T(x, 0) = 0
T_top = 0.5     # T(x, 0.3) = 0.5

# Coordinates of the four interior nodes
coords = {
    'T1': (0.1, 0.2),
    'T2': (0.2, 0.2),
    'T3': (0.1, 0.1),
    'T4': (0.2, 0.1)
}

# Step 4: Formulate the system of linear equations A*T = b
# The general form is:
# 2*α*T(i,j) - λ*T(i,j-1) - λ*T(i,j+1) - T(i-1,j) - T(i+1,j) = -(dx**2)*f(x,y)
# Re-arranging to put knowns on the right side:
# 2*α*T_ij - λ*T_down - λ*T_up - T_left - T_right = -(dx**2)*f
#
# Our unknowns vector is T = [T1, T2, T3, T4]^T
#
# Matrix A coefficients:
# The diagonal element is 2*α.
# The off-diagonal for vertical neighbors (j+1, j-1) is -λ.
# The off-diagonal for horizontal neighbors (i+1, i-1) is -1.

A = np.array([
    # Eq for T1: 4*T1 - 1*T2 - λ*T3 - 0*T4
    [2*alpha_val, -1,          -lambda_val, 0          ],
    # Eq for T2: -1*T1 + 4*T2 - 0*T3 - λ*T4
    [-1,          2*alpha_val, 0,           -lambda_val],
    # Eq for T3: -λ*T1 - 0*T2 + 4*T3 - 1*T4
    [-lambda_val, 0,           2*alpha_val, -1         ],
    # Eq for T4: 0*T1 - λ*T2 - 1*T3 + 4*T4
    [0,           -lambda_val, -1,          2*alpha_val]
])

# Vector b (right-hand side)
# b_k = sum of contributions from boundary neighbors - (dx**2)*f(x_k,y_k)

# RHS for T1 = T(1,2): neighbors are T_left(0,2) and T_top(1,3)
b1 = T_left + T_top - (dx**2) * f(coords['T1'][0], coords['T1'][1])
# RHS for T2 = T(2,2): neighbors are T_right(3,2) and T_top(2,3)
b2 = T_right + T_top - (dx**2) * f(coords['T2'][0], coords['T2'][1])
# RHS for T3 = T(1,1): neighbors are T_left(0,1) and T_bottom(1,0)
b3 = T_left + T_bottom - (dx**2) * f(coords['T3'][0], coords['T3'][1])
# RHS for T4 = T(2,1): neighbors are T_right(3,1) and T_bottom(2,0)
b4 = T_right + T_bottom - (dx**2) * f(coords['T4'][0], coords['T4'][1])

b = np.array([b1, b2, b3, b4])

print("--- Step 3 & 4: System of Equations AT=b ---")
print("Matrix A:")
print(A)
print("\nVector b:")
print(b)

print("\nThe final system of equations is:")
print(f"Eq 1: {A[0,0]:.1f} T1 + ({A[0,1]:.1f}) T2 + ({A[0,2]:.1f}) T3 + {A[0,3]:.1f} T4 = {b[0]:.2f}")
print(f"Eq 2: ({A[1,0]:.1f}) T1 + {A[1,1]:.1f} T2 + {A[1,2]:.1f} T3 + ({A[1,3]:.1f}) T4 = {b[1]:.2f}")
print(f"Eq 3: ({A[2,0]:.1f}) T1 + {A[2,1]:.1f} T2 + {A[2,2]:.1f} T3 + ({A[2,3]:.1f}) T4 = {b[2]:.2f}")
print(f"Eq 4: {A[3,0]:.1f} T1 + ({A[3,1]:.1f}) T2 + ({A[3,2]:.1f}) T3 + {A[3,3]:.1f} T4 = {b[3]:.2f}\n")


# Step 5: Solve the system
T_vector = np.linalg.solve(A, b)

# Step 6: Present the Solution
print("--- Step 5 & 6: Solution ---")
print(f"The solved vector for {{T}} is [T1, T2, T3, T4]:")
print(f"[{T_vector[0]:.4f}, {T_vector[1]:.4f}, {T_vector[2]:.4f}, {T_vector[3]:.4f}]")