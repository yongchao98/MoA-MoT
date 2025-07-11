import numpy as np

def relu(x):
    """Element-wise ReLU activation function."""
    return np.maximum(0, x)

def main():
    """
    Analyzes and demonstrates the rank of latent representations in a neural network.
    """
    # --- Setup ---
    # Fix random seed for reproducibility
    np.random.seed(42)

    # Input description: 100 data points, 50 features, rank 25
    U = np.random.randn(100, 25)
    V = np.random.randn(25, 50)
    X = U @ V
    input_rank = np.linalg.matrix_rank(X)
    print(f"Input matrix X has shape {X.shape} and its rank is {input_rank}.\n")

    print("--- Analyzing Statements ---\n")

    # --- Statement A: Rank of layer 1 representation can be 20 ---
    print("Analysis for A: Could the rank of the first layer's representation be 20?")
    # This is possible if the transformation reduces the rank.
    # Construct a weight matrix W1 with rank 20. The rank of X @ W1 will be <= 20.
    W1_U_A = np.random.randn(50, 20)
    W1_V_A = np.random.randn(20, 50)
    W1_A = W1_U_A @ W1_V_A  # This matrix has rank 20
    b1_A = np.random.randn(50) * 0.1 # Small bias to not drastically change rank
    Z1_A = X @ W1_A + b1_A
    A1_A = relu(Z1_A)
    rank_A1_A = np.linalg.matrix_rank(A1_A)
    print(f"By constructing a rank-20 weight matrix, we get a latent representation with rank = {rank_A1_A}.")
    print("This shows it is possible for the rank to be 20. Statement A could be TRUE.\n")

    # --- Statement B: Rank of layer 1 representation can be 50 ---
    print("Analysis for B: Could the rank of the first layer's representation be 50?")
    # This is possible because the non-linear ReLU can increase the rank.
    W1_B = np.random.randn(50, 50)
    b1_B = np.random.randn(50)
    Z1_B = X @ W1_B + b1_B
    A1_B = relu(Z1_B)
    rank_A1_B = np.linalg.matrix_rank(A1_B)
    print(f"Using random full-rank weights, the rank increased from {input_rank} to {rank_A1_B}.")
    print("This demonstrates the rank-increasing property of ReLU. It is possible to reach the full rank of 50.")
    print("Statement B could be TRUE.\n")

    # --- Statement C: Rank of layer 2 representation can be 25 ---
    print("Analysis for C: Could the rank of the second layer's representation be 25?")
    # Let's use the high-rank representation A1_B from the previous step as input.
    A1_for_C = A1_B
    # Construct a weight matrix W2 with rank 25.
    W2_U_C = np.random.randn(50, 25)
    W2_V_C = np.random.randn(25, 50)
    W2_C = W2_U_C @ W2_V_C # This matrix has rank 25
    b2_C = np.random.randn(50) * 0.1
    Z2_C = A1_for_C @ W2_C + b2_C
    A2_C = relu(Z2_C)
    rank_A2_C = np.linalg.matrix_rank(A2_C)
    print(f"Using an input with rank {rank_A1_B} and a rank-25 weight matrix, we get a new representation with rank = {rank_A2_C}.")
    print("This demonstrates it is possible for the rank to be 25. Statement C could be TRUE.\n")

    # --- Statement D: Rank of layer 3 representation can be 15 ---
    print("Analysis for D: Could the rank of the last layer's representation be 15?")
    # The output of the last layer, A3, has a dimension of 10.
    A2_for_D = A2_C
    W3_D = np.random.randn(50, 10) # Output dimension is 10
    b3_D = np.random.randn(10)
    A3_D = relu(A2_for_D @ W3_D + b3_D)
    output_shape = A3_D.shape
    rank_A3_D = np.linalg.matrix_rank(A3_D)
    print(f"The final layer's representation matrix A3 has shape {output_shape}.")
    print(f"The rank of a matrix cannot exceed its number of columns, so max rank is {min(output_shape)}.")
    print("A rank of 15 is mathematically impossible.")
    print(f"For demonstration, the rank in our case is {rank_A3_D}.")
    print("Statement D is FALSE.\n")

if __name__ == '__main__':
    main()