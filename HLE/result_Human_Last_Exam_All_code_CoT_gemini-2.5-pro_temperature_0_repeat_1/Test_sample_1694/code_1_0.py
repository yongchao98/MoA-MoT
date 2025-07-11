import numpy as np
import cmath

# This script solves the given system of matrix equations to find the sum of the
# first coordinates of the solution matrices X1 and X2.

# --- Step 1: Solve the first equation for X1 ---

# The first equation is:
# A1 * X1^2 + X1^2 * B = C1
# where A1 = [[5, 0], [0, -5]], B = [[6, 0], [0, 6]], C1 = [[-53/12, 0], [0, 0]]
# Since B = 6*I, the equation simplifies to (A1 + 6*I) * X1^2 = C1.

A1 = np.array([[5, 0], [0, -5]])
B = np.array([[6, 0], [0, 6]])
C1 = np.array([[-53/12, 0], [0, 0]])

# Calculate M1 = A1 + B
M1 = A1 + B

# Solve for X1_sq = inv(M1) * C1
try:
    inv_M1 = np.linalg.inv(M1)
    X1_sq = inv_M1 @ C1
except np.linalg.LinAlgError:
    print("Matrix M1 is singular and cannot be inverted.")
    exit()

# The solutions for X1 are of the form [[a, 0], [0, 0]], so a^2 is the (0,0) element of X1_sq.
k1 = X1_sq[0, 0]

# Calculate the two possible values for the first coordinate of X1.
# Since k1 is negative, the roots are imaginary.
sol1_a_1 = cmath.sqrt(k1)
sol1_a_2 = -sol1_a_1

# --- Step 2: Solve the second equation for X2 ---

# The second equation is:
# A2 * X2^2 + X2^2 * B = C2
# where A2 = [[4, 0], [0, -5]], B = [[6, 0], [0, 6]], C2 = [[-3/11, 0], [0, 0]]
# This simplifies to (A2 + 6*I) * X2^2 = C2.

A2 = np.array([[4, 0], [0, -5]])
C2 = np.array([[-3/11, 0], [0, 0]])

# Calculate M2 = A2 + B
M2 = A2 + B

# Solve for X2_sq = inv(M2) * C2
try:
    inv_M2 = np.linalg.inv(M2)
    X2_sq = inv_M2 @ C2
except np.linalg.LinAlgError:
    print("Matrix M2 is singular and cannot be inverted.")
    exit()

# The solutions for X2 are of the form [[a, 0], [0, 0]], so a^2 is the (0,0) element of X2_sq.
k2 = X2_sq[0, 0]

# Calculate the two possible values for the first coordinate of X2.
sol2_a_1 = cmath.sqrt(k2)
sol2_a_2 = -sol2_a_1

# --- Step 3: Calculate and print the total sum ---

# The problem asks for the sum of the first coordinates of the solutions.
# We sum all possible first coordinates found for X1 and X2.
total_sum = sol1_a_1 + sol1_a_2 + sol2_a_1 + sol2_a_2

print("The possible first coordinates for the solution matrix X1 are:")
print(f"x1_a_1 = {sol1_a_1}")
print(f"x1_a_2 = {sol1_a_2}")
print("\nThe possible first coordinates for the solution matrix X2 are:")
print(f"x2_a_1 = {sol2_a_1}")
print(f"x2_a_2 = {sol2_a_2}")

print("\nThe sum of all these first coordinates is calculated as follows:")
print(f"Sum = ({sol1_a_1}) + ({sol1_a_2}) + ({sol2_a_1}) + ({sol2_a_2})")
print(f"Sum = {total_sum}")