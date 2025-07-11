import numpy as np

# 1. Define constants and calculate coefficients
delta_x = 0.1
delta_y = 0.1
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print(f"Calculated Coefficients:")
print(f"λ = (Δx/Δy)² = ({delta_x}/{delta_y})² = {lambda_val:.1f}")
print(f"α = λ + 1 = {lambda_val:.1f} + 1 = {alpha_val:.1f}\n")

# 2. Set up the matrix A and vector B for the system A*T = B
# From the derived equations:
# 4*T1 - 1*T2 - 1*T3 + 0*T4 = 0.48
# -1*T1 + 4*T2 + 0*T3 - 1*T4 = 1.46
# -1*T1 + 0*T2 + 4*T3 - 1*T4 = -0.01
# 0*T1 - 1*T2 - 1*T3 + 4*T4 = 0.98
A = np.array([
    [4, -1, -1,  0],
    [-1, 4,  0, -1],
    [-1, 0,  4, -1],
    [0, -1, -1,  4]
])

B = np.array([0.48, 1.46, -0.01, 0.98])

# 3. Print the system of equations
print("The system of linear equations to solve is A * {T} = {B}:\n")
for i in range(A.shape[0]):
    row_str = " ".join([f"{val:5.1f}" for val in A[i,:]])
    if i == 1:
        print(f"[{row_str}] [T{i+1}] = [{B[i]:6.2f}]")
    else:
        print(f"[{row_str}] [T{i+1}]   [{B[i]:6.2f}]")
print("\n")

# 4. Solve the system for T
try:
    T = np.linalg.solve(A, B)
    # 5. Print the final result
    print("The vector value for {T} is:")
    # The format [T1, T2, T3, T4] is requested
    print(f"[T1, T2, T3, T4] = {T.tolist()}")
    # Final answer in the required format
    print(f"\n<<<{T.tolist()}>>>")

except np.linalg.LinAlgError:
    print("The system of equations could not be solved as the matrix A is singular.")
