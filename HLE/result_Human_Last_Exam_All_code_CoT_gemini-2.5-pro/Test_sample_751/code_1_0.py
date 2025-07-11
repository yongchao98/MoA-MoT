import numpy as np

# Step 1: Define the system matrices
A = np.array([[1, 1, 0],
              [2, 1, 1],
              [0, 2, 0]])

C = np.array([[0, 1, 0],
              [1, 1, 0]])

n = A.shape[0]
print(f"Original A matrix:\n{A}\n")
print(f"Original C matrix:\n{C}\n")

# Step 2: Construct the transformation matrix T
# We need to find n=3 linearly independent rows from C, CA, ...
# Let c1 and c2 be the rows of C
c1 = C[0, :]
c2 = C[1, :]

# Calculate rows of CA
c1A = c1 @ A
c2A = c2 @ A

# Select n=3 linearly independent rows to form the transformation matrix T.
# The standard procedure checks for linear independence in the order c1, c2, c1A, c2A...
# 1. c1 = [0, 1, 0] is independent.
# 2. c2 = [1, 1, 0] is linearly independent from c1.
# 3. c1A = [2, 1, 1] is linearly independent from c1 and c2.
# We have found 3 linearly independent rows.
# The convention is to group rows related to the same output together in T.
# Since we used c1 and c1A (related to the first output) and c2 (related to the second),
# the observability indices are mu_1 = 2, mu_2 = 1.
# We form T by stacking these rows:
T = np.vstack([c1, c1A, c2])

print(f"Transformation matrix T:\n{T}\n")

# Step 3: Calculate the inverse of T
try:
    T_inv = np.linalg.inv(T)
    print(f"Inverse of T (T_inv):\n{T_inv}\n")
except np.linalg.LinAlgError:
    print("Matrix T is singular, cannot compute inverse.")
    exit()

# Step 4: Compute the new output matrix C_o for the observer canonical form
# C_o = C * T_inv
C_o = C @ T_inv

# Clean up the result to show integer values, as they are exact in this case
C_o_clean = np.round(C_o).astype(int)

print("The new matrix C in observer canonical form (C_o) is:")
# Print the final matrix element by element as requested
print("[", end="")
for i, row in enumerate(C_o_clean):
    if i > 0:
        print(" ", end="")
    print("[", end="")
    for j, val in enumerate(row):
        print(f"{val}", end="")
        if j < len(row) - 1:
            print(", ", end="")
    print("]", end="")
    if i < len(C_o_clean) - 1:
        print(",\n", end="")
print("]")

# The final answer in the required format
final_answer = C_o_clean.tolist()
# <<<[[1, 0, 0], [0, 0, 1]]>>>