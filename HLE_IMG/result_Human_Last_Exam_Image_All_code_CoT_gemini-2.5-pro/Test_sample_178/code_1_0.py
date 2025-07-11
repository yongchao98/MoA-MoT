import numpy as np

# Step 1: Define grid parameters and the source function
# The domain is a 0.3 x 0.3 square with 3 intervals in each direction.
dx = 0.3 / 3.0
dy = 0.3 / 3.0

# Source function f(x, y) = 100xy
def f(x, y):
    """Calculates the value of the source function f(x,y)."""
    return 100 * x * y

# Boundary conditions are given as constants
T_left = 0.0   # T(0, y) = 0
T_right = 1.0  # T(0.3, y) = 1
T_bottom = 0.0 # T(x, 0) = 0
T_top = 0.5    # T(x, 0.3) = 0.5

# Step 2: Calculate coefficients lambda and alpha
lambda_val = (dx / dy)**2
alpha = lambda_val + 1

print("1. Calculated Coefficients:")
print(f"   Δx = {dx}")
print(f"   Δy = {dy}")
print(f"   λ = (Δx/Δy)² = {lambda_val:.1f}")
print(f"   α = λ + 1 = {alpha:.1f}\n")

# The finite difference equation is:
# 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)² * f(x,y)
# With λ=1 and α=2, it simplifies to:
# 4*T(i,j) - T(i,j+1) - T(i,j-1) - T(i+1,j) - T(i-1,j) = -(Δx)² * f(x,y)

# Step 3: Set up the system of linear equations A * {T} = {b}
# The unknowns T1, T2, T3, T4 are at coordinates:
# T1: (x=0.1, y=0.2), T2: (x=0.2, y=0.2)
# T3: (x=0.1, y=0.1), T4: (x=0.2, y=0.1)
# The vector of unknowns is {T} = [T1, T2, T3, T4]

# The coefficient matrix A is derived from the left side of the equations.
A = np.array([
    # Eq for T1: 4*T1 - T2 - T3 = ...
    [4, -1, -1,  0],
    # Eq for T2: 4*T2 - T1 - T4 = ...
    [-1, 4,  0, -1],
    # Eq for T3: 4*T3 - T1 - T4 = ...
    [-1, 0,  4, -1],
    # Eq for T4: 4*T4 - T2 - T3 = ...
    [0, -1, -1,  4]
])

# The constant vector b is derived by moving known terms to the right side.
# b_k = sum_of_known_neighbor_temps - (dx**2) * f(x,y)

# RHS for T1 equation
b1 = T_top + T_left - (dx**2) * f(0.1, 0.2)

# RHS for T2 equation
b2 = T_top + T_right - (dx**2) * f(0.2, 0.2)

# RHS for T3 equation
b3 = T_bottom + T_left - (dx**2) * f(0.1, 0.1)

# RHS for T4 equation
b4 = T_bottom + T_right - (dx**2) * f(0.2, 0.1)

b = np.array([b1, b2, b3, b4])

print("2. System of Linear Equations (A * {T} = {b}):")
# This loop prints each equation with its numerical coefficients.
labels = ["T1", "T2", "T3", "T4"]
for i in range(4):
    eq_str = ""
    for j in range(4):
        eq_str += f"{A[i, j]:+5.2f}*{labels[j]} "
    eq_str += f"= {b[i]:.4f}"
    print(f"   {eq_str}")
print("")

# Step 4: Solve the system of equations for {T}
T_solution = np.linalg.solve(A, b)

# Step 5: Print the final result
print("3. Final Solution Vector {T}:")
print(f"   [T1, T2, T3, T4] = [{T_solution[0]:.6f}, {T_solution[1]:.6f}, {T_solution[2]:.6f}, {T_solution[3]:.6f}]")

final_answer_str = f"[{T_solution[0]:.6f}, {T_solution[1]:.6f}, {T_solution[2]:.6f}, {T_solution[3]:.6f}]"
# print(f"\n<<<{final_answer_str}>>>")