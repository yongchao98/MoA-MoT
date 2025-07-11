import numpy as np

# Step 1 & 2: Define parameters and calculate coefficients
delta_x = 0.3 / 3.0
delta_y = 0.3 / 3.0
lmbda = (delta_x / delta_y)**2
alpha = lmbda + 1

print(f"Calculated coefficients:")
print(f"lambda = ({delta_x}/{delta_y})^2 = {lmbda}")
print(f"alpha = lambda + 1 = {alpha}\n")

# Step 3: Formulate the system of linear equations A * {T} = {B}
# The general equation becomes:
# 4*T(i,j) - T(i,j+1) - T(i,j-1) - T(i+1,j) - T(i-1,j) = -x(i)*y(j)

# Boundary values
T_top = 0.5     # T(x, 0.3)
T_bottom = 0    # T(x, 0)
T_left = 0      # T(0, y)
T_right = 1.0   # T(0.3, y)

# Node coordinates
x = [0, 0.1, 0.2, 0.3]
y = [0, 0.1, 0.2, 0.3]

# Equation for T1 at (x1, y2) = (0.1, 0.2)
# 4*T1 - T(x1,y3) - T3 - T2 - T(x0,y2) = -x1*y2
# 4*T1 - T_top - T3 - T2 - T_left = -0.1*0.2
# 4*T1 - T2 - T3 = T_top + T_left - 0.02 = 0.5 + 0 - 0.02 = 0.48
b1 = 0.48

# Equation for T2 at (x2, y2) = (0.2, 0.2)
# 4*T2 - T(x2,y3) - T4 - T(x3,y2) - T1 = -x2*y2
# 4*T2 - T_top - T4 - T_right - T1 = -0.2*0.2
# -T1 + 4*T2 - T4 = T_top + T_right - 0.04 = 0.5 + 1.0 - 0.04 = 1.46
b2 = 1.46

# Equation for T3 at (x1, y1) = (0.1, 0.1)
# 4*T3 - T1 - T(x1,y0) - T4 - T(x0,y1) = -x1*y1
# 4*T3 - T1 - T_bottom - T4 - T_left = -0.1*0.1
# -T1 + 4*T3 - T4 = T_bottom + T_left - 0.01 = 0 + 0 - 0.01 = -0.01
b3 = -0.01

# Equation for T4 at (x2, y1) = (0.2, 0.1)
# 4*T4 - T2 - T(x2,y0) - T(x3,y1) - T3 = -x2*y1
# 4*T4 - T2 - T_bottom - T_right - T3 = -0.2*0.1
# -T2 - T3 + 4*T4 = T_bottom + T_right - 0.02 = 0 + 1.0 - 0.02 = 0.98
b4 = 0.98

# Coefficient Matrix A
A = np.array([
    [ 4, -1, -1,  0], # Equation for T1
    [-1,  4,  0, -1], # Equation for T2
    [-1,  0,  4, -1], # Equation for T3
    [ 0, -1, -1,  4]  # Equation for T4
])

# Constant vector B
B = np.array([b1, b2, b3, b4])

print("The system of linear equations A*{T} = {B} is:")
print("where {T} = [T1, T2, T3, T4]^T\n")
print("Coefficient matrix A:")
print(A)
print("\nConstant vector B:")
print(B.reshape(-1, 1))

# Step 4: Solve the system
T_vector = np.linalg.solve(A, B)

# Print the final result
print("\nThe solution vector for {T} is:")
print(f"[{T_vector[0]}, {T_vector[1]}, {T_vector[2]}, {T_vector[3]}]")
<<<[0.3016666666666667, 0.5470833333333333, 0.17958333333333333, 0.4266666666666667]>>>