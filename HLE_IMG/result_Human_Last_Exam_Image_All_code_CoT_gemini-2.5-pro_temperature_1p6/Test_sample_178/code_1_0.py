import numpy as np

# 1. Define grid parameters and calculate coefficients
num_intervals = 3
L = 0.3
dx = L / num_intervals
dy = L / num_intervals
lambda_val = (dx / dy)**2
alpha_val = lambda_val + 1

print("--- Step 1: Computed Coefficients ---")
print(f"λ = (Δx/Δy)² = ({dx}/{dy})² = {lambda_val:.2f}")
print(f"α = λ + 1 = {lambda_val:.2f} + 1 = {alpha_val:.2f}")
print("-" * 35)

# 2. Define the system of linear equations AT = b derived from the finite difference method
# The unknowns are T = [T1, T2, T3, T4]^T

# Equation for T1: 4*T1 - T2 - T3 = 0.48
# Equation for T2: -T1 + 4*T2 - T4 = 1.46
# Equation for T3: -T1 + 4*T3 - T4 = -0.01
# Equation for T4: -T2 - T3 + 4*T4 = 0.98

print("\n--- Step 2: System of Linear Equations (AT = b) ---")
# Coefficient matrix A
A = np.array([
    [ 4, -1, -1,  0],  # Equation for T1
    [-1,  4,  0, -1],  # Equation for T2
    [-1,  0,  4, -1],  # Equation for T3
    [ 0, -1, -1,  4]   # Equation for T4
])

# Constant vector b
b = np.array([
    0.48,   # Right-hand side for T1 equation
    1.46,   # Right-hand side for T2 equation
    -0.01,  # Right-hand side for T3 equation
    0.98    # Right-hand side for T4 equation
])

print("Matrix A:")
print(A)
print("\nVector b:")
print(b.reshape(-1, 1))
print("-" * 35)


# 3. Solve the system for T
T = np.linalg.solve(A, b)

print("\n--- Step 3: Final Solution ---")
print("The vector value for {T} is:")
print(f"[T1, T2, T3, T4] = [{T[0]:.6f}, {T[1]:.6f}, {T[2]:.6f}, {T[3]:.6f}]")

# Final Answer Block
answer = f"[{T[0]}, {T[1]}, {T[2]}, {T[3]}]"
# print(f"\n<<<{answer}>>>")