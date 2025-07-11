import numpy as np

# Step 1: Define physical and numerical parameters
delta_x = 0.3 / 3.0
delta_y = 0.3 / 3.0

# Step 2: Calculate coefficients lambda and alpha
lambda_val = (delta_x / delta_y)**2
alpha_val = lambda_val + 1

print("1. Calculated Coefficients:")
print(f"Δx = {delta_x}")
print(f"Δy = {delta_y}")
print(f"λ = (Δx/Δy)² = ({delta_x}/{delta_y})² = {lambda_val:.2f}")
print(f"α = λ + 1 = {lambda_val:.2f} + 1 = {alpha_val:.2f}")
print("-" * 40)

# Step 3 & 4: Formulate the system of equations by applying the finite difference
# equation to each node and incorporating boundary conditions.
# The general equation is:
# 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)²*f(x,y)
# With λ=1, α=2, f(x,y)=100xy, and Δx=0.1, this simplifies to:
# 4*T_center - T_up - T_down - T_right - T_left = -x*y

# Node T1(x=0.1, y=0.2): 4*T1 - T(x=0.1,y=0.3) - T3 - T2 - T(x=0,y=0.2) = -(0.1*0.2)
# 4*T1 - 0.5 - T3 - T2 - 0 = -0.02  =>  4*T1 - T2 - T3 = 0.48

# Node T2(x=0.2, y=0.2): 4*T2 - T(x=0.2,y=0.3) - T4 - T(x=0.3,y=0.2) - T1 = -(0.2*0.2)
# 4*T2 - 0.5 - T4 - 1.0 - T1 = -0.04  =>  -T1 + 4*T2 - T4 = 1.46

# Node T3(x=0.1, y=0.1): 4*T3 - T1 - T(x=0.1,y=0) - T4 - T(x=0,y=0.1) = -(0.1*0.1)
# 4*T3 - T1 - 0 - T4 - 0 = -0.01  =>  -T1 + 4*T3 - T4 = -0.01

# Node T4(x=0.2, y=0.1): 4*T4 - T2 - T(x=0.2,y=0) - T(x=0.3,y=0.1) - T3 = -(0.2*0.1)
# 4*T4 - T2 - 0 - 1.0 - T3 = -0.02  =>  -T2 - T3 + 4*T4 = 0.98

# Step 5: Set up the matrix A and vector B for the system A*T = B
A = np.array([
    [4.0, -1.0, -1.0,  0.0],
    [-1.0, 4.0,  0.0, -1.0],
    [-1.0, 0.0,  4.0, -1.0],
    [0.0, -1.0, -1.0,  4.0]
])

B = np.array([0.48, 1.46, -0.01, 0.98])

print("2. The system of linear equations (A*T = B) is:")
print(f"Eq 1: {A[0,0]:.1f}*T1 + ({A[0,1]:.1f})*T2 + ({A[0,2]:.1f})*T3 + {A[0,3]:.1f}*T4 = {B[0]}")
print(f"Eq 2: ({A[1,0]:.1f})*T1 + {A[1,1]:.1f}*T2 + {A[1,2]:.1f}*T3 + ({A[1,3]:.1f})*T4 = {B[1]}")
print(f"Eq 3: ({A[2,0]:.1f})*T1 + {A[2,1]:.1f}*T2 + {A[2,2]:.1f}*T3 + ({A[2,3]:.1f})*T4 = {B[2]}")
print(f"Eq 4: {A[3,0]:.1f}*T1 + ({A[3,1]:.1f})*T2 + ({A[3,2]:.1f})*T3 + {A[3,3]:.1f}*T4 = {B[3]}")
print("-" * 40)

# Step 6: Solve the system for T
T = np.linalg.solve(A, B)

# Step 7: Print the final result
print("3. The resulting vector value for {T} is:")
print(f"T1 = {T[0]:.4f}")
print(f"T2 = {T[1]:.4f}")
print(f"T3 = {T[2]:.4f}")
print(f"T4 = {T[3]:.4f}")
print("\nIn [T1, T2, T3, T4] format:")
print(f"[{T[0]:.4f}, {T[1]:.4f}, {T[2]:.4f}, {T[3]:.4f}]")