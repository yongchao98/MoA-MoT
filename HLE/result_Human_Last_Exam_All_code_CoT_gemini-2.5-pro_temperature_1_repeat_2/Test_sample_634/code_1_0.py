import numpy as np

def print_poly_matrix(matrix):
    """Helper function to print a matrix of polynomials nicely."""
    max_lens = [0] * matrix.shape[1]
    for j in range(matrix.shape[1]):
        for i in range(matrix.shape[0]):
            max_lens[j] = max(max_lens[j], len(str(matrix[i, j]).replace('\n', '')))

    print("[", end="")
    for i in range(matrix.shape[0]):
        if i > 0:
            print(" ", end="")
        print("[", end="")
        for j in range(matrix.shape[1]):
            s = str(matrix[i, j]).replace('\n', '')
            print(s.ljust(max_lens[j]), end="")
            if j < matrix.shape[1] - 1:
                print(", ", end="")
        print("]", end="")
        if i < matrix.shape[0] - 1:
            print(",\n")
    print("]")


# Use numpy.poly1d for polynomial arithmetic and representation
s = np.poly1d([1, 0])

# Define the polynomial matrices P1 and P2
P1 = np.array([
    [s**2 + s, -s],
    [-s**2 - np.poly1d(1), s**2]
], dtype=object)

P2 = np.array([
    [s, np.poly1d(0)],
    [-s - np.poly1d(1), np.poly1d(1)]
], dtype=object)

# Form the 4x2 stacked matrix M
M = np.vstack([P1, P2])

# Perform elementary row operations to reduce the matrix M.
# Let the rows be R1, R2, R3, R4.

# Step 1: Use R3 to simplify R1, R2, and R4.
# R1_new = R1 - (s+1)*R3
M[0] = M[0] - (s + 1) * M[2]
# R2_new = R2 + s*R3
M[1] = M[1] + s * M[2]
# R4_new = R4 + R3
M[3] = M[3] + M[2]

# Step 2: Use the new R4 to simplify the new R2 and R3.
# R2_new = R2 - R4
M[1] = M[1] - M[3]
# R3_new = R3 + s*R4
M[2] = M[2] + s * M[3]

# Step 3: Use the new R3 to simplify the new R1 and R2.
# R1_new = R1 + R3
M[0] = M[0] + M[2]
# R2_new = R2 - s*R3
M[1] = M[1] - s * M[2]

# Step 4: Use the new R2 to simplify the new R3 and R4.
# R3_new = R3 + s*R2
M[2] = M[2] + s * M[1]
# R4_new = R4 + R2
M[3] = M[3] + M[1]

# After reduction, the non-zero rows form the GCRD.
# Extract the non-zero rows to form the GCRD matrix.
non_zero_rows = []
for row in M:
    if not all(p == np.poly1d(0) for p in row):
        non_zero_rows.append(row)

# The result of the row reduction is [[0,-1], [-1,0]].
# We arrange them to be in a lower-triangular form.
gcrd_matrix = np.array([non_zero_rows[1], non_zero_rows[0]], dtype=object)

print("The calculated greatest common right divisor is:")
print_poly_matrix(gcrd_matrix)

print("\n\nThis matrix is unimodular (determinant is a non-zero constant),")
print("and is equivalent to the identity matrix, which is the canonical GCRD.")
print("This implies that the original matrices P1 and P2 are right co-prime.")
print("The final canonical GCRD is:")

I = np.array([[np.poly1d([1]), np.poly1d([0])],
              [np.poly1d([0]), np.poly1d([1])]])

print_poly_matrix(I)
