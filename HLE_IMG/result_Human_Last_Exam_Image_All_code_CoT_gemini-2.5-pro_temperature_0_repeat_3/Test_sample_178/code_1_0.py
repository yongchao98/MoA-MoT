import numpy as np

# Step 1 & 2: Define parameters
# The problem is discretized into a grid. The unknowns are T1, T2, T3, T4.
# We will set up a system of 4 linear equations A*T = B, where T = [T1, T2, T3, T4].

# Step 3 & 4: Formulate the system of equations
# The general equation is: 4*T_ij - T_i,j+1 - T_i,j-1 - T_i+1,j - T_i-1,j = -x_i * y_j

# Equation for T1 (x=0.1, y=0.2)
# 4*T1 - T(x=0.1,y=0.3) - T3 - T2 - T(x=0,y=0.2) = -0.1 * 0.2
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.02
# 4*T1 - T2 - T3 = 0.48

# Equation for T2 (x=0.2, y=0.2)
# 4*T2 - T(x=0.2,y=0.3) - T4 - T(x=0.3,y=0.2) - T1 = -0.2 * 0.2
# 4*T2 - 0.5 - T4 - 1 - T1 = -0.04
# -T1 + 4*T2 - T4 = 1.46

# Equation for T3 (x=0.1, y=0.1)
# 4*T3 - T1 - T(x=0.1,y=0) - T4 - T(x=0,y=0.1) = -0.1 * 0.1
# 4*T3 - T1 - 0 - T4 - 0 = -0.01
# -T1 + 4*T3 - T4 = -0.01

# Equation for T4 (x=0.2, y=0.1)
# 4*T4 - T2 - T(x=0.2,y=0) - T(x=0.3,y=0.1) - T3 = -0.2 * 0.1
# 4*T4 - T2 - 0 - 1 - T3 = -0.02
# -T2 - T3 + 4*T4 = 0.98

# Step 5: Solve the system using NumPy
# The system in matrix form A * {T} = {B} is:
# T = [T1, T2, T3, T4]
A = np.array([
    [ 4, -1, -1,  0],
    [-1,  4,  0, -1],
    [-1,  0,  4, -1],
    [ 0, -1, -1,  4]
])

B = np.array([0.48, 1.46, -0.01, 0.98])

# Print the system of equations
print("The system of linear equations is:")
for i in range(A.shape[0]):
    equation = []
    for j in range(A.shape[1]):
        if A[i, j] != 0:
            equation.append(f"{A[i, j]:.0f} * T{j+1}")
    print(" + ".join(equation).replace("+ -", "- "), f"= {B[i]}")

# Solve for T
T = np.linalg.solve(A, B)

# Print the final result
print("\nThe vector value for {T} is:")
print(f"[{T[0]}, {T[1]}, {T[2]}, {T[3]}]")

# Final answer in the required format
answer = f"[{T[0]}, {T[1]}, {T[2]}, {T[3]}]"
# print(f"\n<<<{answer}>>>")