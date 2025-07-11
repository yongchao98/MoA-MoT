import numpy as np

# This script solves the system of linear equations derived from the
# finite difference approximation of the given Poisson equation.

# 1. Define grid parameters
L = 0.3  # Side length of the square plate
intervals = 3
dx = L / intervals
dy = L / intervals

# 2. Calculate coefficients
lambda_val = (dx / dy)**2
alpha = lambda_val + 1

print("Calculated Coefficients:")
print(f"λ = {lambda_val:.1f}")
print(f"α = {alpha:.1f}")
print("-" * 50)

# 3. Formulate the system of equations A * T = B
# The four equations for T1, T2, T3, T4 are derived by applying the
# finite difference formula at each interior node and substituting the boundary conditions.
# The general form is: 4*T_ij - T_neighbors = -x*y

A = np.array([
    [4, -1, -1,  0],  # Equation for T1
    [-1, 4,  0, -1],  # Equation for T2
    [-1, 0,  4, -1],  # Equation for T3
    [0, -1, -1,  4]   # Equation for T4
])

B = np.array([
    0.48,   # RHS for T1: 0.5 (T_top) + 0.0 (T_left) - (-0.1*0.2)
    1.46,   # RHS for T2: 0.5 (T_top) + 1.0 (T_right) - (-0.2*0.2)
    -0.01,  # RHS for T3: 0.0 (T_bottom) + 0.0 (T_left) - (-0.1*0.1)
    0.98    # RHS for T4: 0.0 (T_bottom) + 1.0 (T_right) - (-0.2*0.1)
])

# Print the final equations with all numerical coefficients
print("The final system of linear equations is:")
variables = ['T1', 'T2', 'T3', 'T4']
for i in range(len(B)):
    parts = []
    for j in range(len(A[i])):
        coeff = A[i, j]
        # Only include terms with non-zero coefficients
        if coeff != 0:
            parts.append(f"({int(coeff)})*{variables[j]}")
    
    # Join the parts with '+' and clean up for negative signs
    equation_str = " + ".join(parts).replace("+ (-", "- (")
    print(f"Eq {i+1}: {equation_str} = {B[i]}")

print("-" * 50)

# 4. Solve the system for the temperature vector T
try:
    T = np.linalg.solve(A, B)
    # Print the final result in the specified format
    print("The solution vector for {T} is:")
    print(f"[T1, T2, T3, T4] = [{T[0]:.5f}, {T[1]:.5f}, {T[2]:.5f}, {T[3]:.5f}]")
except np.linalg.LinAlgError:
    print("The system of equations could not be solved.")
