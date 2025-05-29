import numpy as np

# Initial matrix
matrix = np.array([
    [1, 8, 1, 5, 8],
    [4, 4, 5, 7, 0],
    [8, 8, 2, 2, 0],
    [4, 2, 7, 9, 3],
    [2, 2, 2, 8, 3],
    [2, 6, 0, 3, 8]
])

# Step 2: Set all elements divisible by 5 to zero
matrix[matrix % 5 == 0] = 0

# Step 3: Map each occurrence of 3 to 4
matrix[matrix == 3] = 4

# Step 4: Remove every 5-th row (1-indexed)
matrix = np.delete(matrix, 4, axis=0)

# Step 5: Remove every 3-th row (1-indexed)
matrix = np.delete(matrix, 2, axis=0)

# Step 6: Mirror the matrix along the counterdiagonal
n = matrix.shape[0]
m = matrix.shape[1]
mirrored_matrix = np.zeros((m, n), dtype=int)
for i in range(n):
    for j in range(m):
        mirrored_matrix[j, i] = matrix[n - i - 1, m - j - 1]

# Step 7: Set all elements divisible by 3 to zero
mirrored_matrix[mirrored_matrix % 3 == 0] = 0

# Print the final matrix
print(mirrored_matrix)