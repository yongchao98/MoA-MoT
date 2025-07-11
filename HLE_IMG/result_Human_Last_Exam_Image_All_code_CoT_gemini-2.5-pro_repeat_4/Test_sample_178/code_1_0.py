import numpy as np

# Step 1 & 2: Define parameters and calculate coefficients
delta_x = 0.1
delta_y = 0.1
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print(f"Calculated lambda (λ) = {lambda_val}")
print(f"Calculated alpha (α) = {alpha_val}\n")

# Step 3 & 4: Formulate the system of linear equations A*T = B

# The four equations are:
# For T1 (x=0.1, y=0.2): 4*T1 - T2 - T3 = 0.5 - (0.1*0.2) = 0.48
# For T2 (x=0.2, y=0.2): 4*T2 - T1 - T4 = 1 + 0.5 - (0.2*0.2) = 1.46
# For T3 (x=0.1, y=0.1): 4*T3 - T1 - T4 = 0 + 0 - (0.1*0.1) = -0.01
# For T4 (x=0.2, y=0.1): 4*T4 - T2 - T3 = 1 + 0 - (0.2*0.1) = 0.98

# Let T = [T1, T2, T3, T4]
# The system is:
# 4*T1 - 1*T2 - 1*T3 + 0*T4 = 0.48
# -1*T1 + 4*T2 + 0*T3 - 1*T4 = 1.46
# -1*T1 + 0*T2 + 4*T3 - 1*T4 = -0.01
# 0*T1 - 1*T2 - 1*T3 + 4*T4 = 0.98

# Define the coefficient matrix A
A = np.array([
    [4, -1, -1, 0],
    [-1, 4, 0, -1],
    [-1, 0, 4, -1],
    [0, -1, -1, 4]
])

# Define the constant vector B
B = np.array([0.48, 1.46, -0.01, 0.98])

# Print the system of equations
print("The system of linear equations is A * {T} = {B}:\n")
variables = ['T1', 'T2', 'T3', 'T4']
for i in range(A.shape[0]):
    equation = ""
    for j in range(A.shape[1]):
        if A[i, j] != 0:
            equation += f"({A[i, j]}){variables[j]} + "
    equation = equation[:-3]  # Remove trailing ' + '
    equation += f" = {B[i]}"
    print(f"Equation {i+1}: {equation}")

# Step 5: Solve the system for T
T_solution = np.linalg.solve(A, B)

# Print the final result
print("\n----------------------------------------------------")
print(f"The vector value for {{T}} is: {T_solution.tolist()}")
print("----------------------------------------------------")

# Final answer in the requested format
final_answer = T_solution.tolist()
print(f'<<<{final_answer}>>>')