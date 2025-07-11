import numpy as np

# Step 1: Find the roots of the characteristic equation 8r^2 - 6r + 1 = 0.
# The roots can be found using the quadratic formula or by factoring: (4r - 1)(2r - 1) = 0.
# The roots are r1 = 1/2 and r2 = 1/4.
# We will assign the larger root to B and the smaller one to D.
B = 0.5
D = 0.25

# Step 2: Find the particular solution (E).
# Assume a constant solution y_p[n] = K for the constant input of 1.
# 8K - 6K + K = 1  =>  3K = 1  =>  K = 1/3.
# So, E is 1/3.
E = 1/3.0

# Step 3: Use initial conditions to set up equations for A and C.
# The general solution is y[n] = A*B^n + C*D^n + E.
# From y[0] = 1:  A*(B^0) + C*(D^0) + E = 1  => A + C + 1/3 = 1 => A + C = 2/3.
# From y[-1] = 2: A*(B^-1) + C*(D^-1) + E = 2 => A/B + C/D + E = 2 => 2A + 4C + 1/3 = 2 => 2A + 4C = 5/3.

# Step 4: Solve the system of linear equations for A and C.
# Eq1: A + C = 2/3
# Eq2: 2A + 4C = 5/3
# We can represent this system as M*x = v, where x = [A, C]
M = np.array([[1, 1], [2, 4]])
v = np.array([2/3.0, 5/3.0])
solution = np.linalg.solve(M, v)
A = solution[0]
C = solution[1]

# The final equation form is y[n] = A * B^n + C * D^n + E
print(f"The parameters for the final equation are:")
print(f"A = {A} (or 1/2)")
print(f"B = {B} (or 1/2)")
print(f"C = {C} (or 1/6)")
print(f"D = {D} (or 1/4)")
print(f"E = {E} (or 1/3)")
print("")

# Step 5: Calculate the final expression: E/A + (D*C)/B
result = (E / A) + (D * C) / B

print("The result of the expression E/A + (D*C)/B is:")
print(result)