import numpy as np
from fractions import Fraction

# Step 1: Find roots of the characteristic equation 8r^2 - 6r + 1 = 0
a, b, c = 8, -6, 1
delta = (b**2 - 4*a*c)**0.5
r1 = (-b + delta) / (2*a)
r2 = (-b - delta) / (2*a)

# Step 2: Find the particular solution constant E
# 8E - 6E + E = 1 => 3E = 1
E = 1/3

# Step 3: Assign B and D, assuming B is the root with the larger magnitude
if abs(r1) > abs(r2):
    B = r1
    D = r2
else:
    B = r2
    D = r1

# Step 4: Use initial conditions to create a system of equations for A and C
# y[0] = 1 => A + C + E = 1
# y[-1] = 2 => A/B + C/D + E = 2
y0 = 1
y_minus_1 = 2

# M * [A, C]^T = V
M = np.array([
    [1, 1],
    [1/B, 1/D]
])
V = np.array([
    y0 - E,
    y_minus_1 - E
])

# Step 5: Solve for A and C
solution = np.linalg.solve(M, V)
A = solution[0]
C = solution[1]

# Step 6: Calculate the final expression
result = E/A + (D*C)/B

# Display the components of the expression and the final result
# The requested expression is E/A + (D*C)/B
print("The values of the parameters are:")
print(f"E = {E:.4f} ({Fraction(E).limit_denominator()})")
print(f"A = {A:.4f} ({Fraction(A).limit_denominator()})")
print(f"D = {D:.4f} ({Fraction(D).limit_denominator()})")
print(f"C = {C:.4f} ({Fraction(C).limit_denominator()})")
print(f"B = {B:.4f} ({Fraction(B).limit_denominator()})")
print("\nCalculating the expression E/A + (D*C)/B:")
print(f"({E}) / ({A}) + (({D}) * ({C})) / ({B}) = {result}")
