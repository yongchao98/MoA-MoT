import numpy as np

def relu(x):
    """Element-wise ReLU activation function."""
    return np.maximum(0, x)

def get_rank(matrix):
    """Computes the rank of a matrix."""
    return np.linalg.matrix_rank(matrix)

# For reproducibility
np.random.seed(42)

# --- Setup: Input Matrix ---
# Description: 100 data points, 50 features, rank 25
num_datapoints = 100
input_features = 50
input_rank = 25

# Create an input matrix X with the specified shape and rank
U = np.random.randn(num_datapoints, input_rank)
V = np.random.randn(input_rank, input_features)
X = U @ V
print(f"--- Input Data ---")
print(f"Input matrix X has shape {X.shape} and rank {get_rank(X)} (expected {input_rank}).\n")


# --- Network Layer Dimensions ---
dim1 = 50
dim2 = 50
dim3 = 10

# --- Analysis of Statement A: rank(A1) = 20 (Plausible) ---
print("--- Analysis of Statement A: rank(A1) = 20 ---")
# A large negative bias can zero out many activations, causing the rank to decrease.
# We can construct weights/biases to achieve this.
W1_A = np.random.randn(input_features, dim1)
b1_A = np.ones(dim1) * -2.5  # A carefully chosen negative bias
Z1_A = X @ W1_A + b1_A
A1_A = relu(Z1_A)
rank_A1_A = get_rank(A1_A)
print(f"By causing some activations to become zero, rank can decrease.")
print(f"For one such choice of W1 and b1, rank(A1) = {rank_A1_A}.")
print("It is plausible to find parameters that result in a rank of exactly 20.")
print("Statement A could be True.\n")


# --- Analysis of Statement B: rank(A1) = 50 (Plausible) ---
print("--- Analysis of Statement B: rank(A1) = 50 ---")
# The ReLU non-linearity can increase rank, up to the maximum possible for the matrix shape.
W1_B = np.random.randn(input_features, dim1)
b1_B = np.random.randn(dim1) * 0.1  # Small random bias
Z1_B = X @ W1_B + b1_B
A1_B = relu(Z1_B)
rank_A1_B = get_rank(A1_B)
rank_Z1_B = get_rank(Z1_B)
print(f"The rank before ReLU, rank(Z1), is {rank_Z1_B} (at most {input_rank}+1 = {input_rank + 1}).")
print(f"After ReLU, the non-linearity can increase the rank.")
print(f"With random W1 and b1, rank(A1) = {rank_A1_B}.")
print(f"The rank can reach the maximum possible value of min({num_datapoints}, {dim1}) = {min(num_datapoints, dim1)}.")
print("Statement B could be True.\n")


# --- Analysis of Statement C: rank(A2) = 25 (Plausible) ---
print("--- Analysis of Statement C: rank(A2) = 25 ---")
# We can construct a weight matrix W2 with rank 25 to constrain the output rank.
# We use the high-rank A1_B from the previous step as input.
A1_input = A1_B
# Create a weight matrix W2 with rank 25
W2_C = np.random.randn(dim1, 25) @ np.random.randn(25, dim2)
b2_C = np.random.randn(dim2)
Z2_C = A1_input @ W2_C + b2_C
A2_C = relu(Z2_C)
rank_A2_C = get_rank(A2_C)
print(f"Using an input A1 with rank {get_rank(A1_input)} and a W2 with rank {get_rank(W2_C)},")
print(f"the resulting rank of the second latent matrix A2 is {rank_A2_C}.")
print("It is plausible for the rank to be exactly 25. Statement C could be True.\n")


# --- Analysis of Statement D: rank(A3) = 15 (Impossible) ---
print("--- Analysis of Statement D: rank(A3) = 15 ---")
# The final latent matrix A3 has 10 columns. Its rank cannot exceed 10.
A2_input = A2_C
W3_D = np.random.randn(dim2, dim3)
b3_D = np.random.randn(dim3)
A3_D = relu(A2_input @ W3_D + b3_D)
rank_A3_D = get_rank(A3_D)
print(f"The final latent matrix A3 has shape {A3_D.shape}.")
print(f"The rank of a matrix is at most the minimum of its dimensions.")
print(f"max_rank(A3) = min({A3_D.shape[0]}, {A3_D.shape[1]}) = {min(A3_D.shape)}.")
print(f"A rank of 15 is therefore impossible. A sample calculation gives rank(A3) = {rank_A3_D}.")
print("Statement D cannot be True.\n")
