import numpy as np

def solve():
    """
    Analyzes and simulates the rank of latent representations in a 3-layer MLP.
    """
    # --- Setup ---
    # Network dimensions and input properties
    N_DATAPOINTS = 100
    D_FEATURES = 50
    D_HIDDEN = 50
    D_OUTPUT = 10
    INPUT_RANK = 25

    # Set a seed for reproducibility of the simulation
    np.random.seed(42)

    # --- 1. Create Input Data ---
    # Create an input matrix X with shape (100, 50) and rank 25.
    # We do this by multiplying two smaller random matrices.
    # A: (100, 25), B: (25, 50). The product X = A @ B will have rank at most 25.
    # With random values, the rank will almost certainly be exactly 25.
    A = np.random.randn(N_DATAPOINTS, INPUT_RANK)
    B = np.random.randn(INPUT_RANK, D_FEATURES)
    X = A @ B
    rank_X = np.linalg.matrix_rank(X)

    # --- 2. Define Network Layers ---
    # We assume weights are initialized randomly and biases are non-zero.
    W1 = np.random.randn(D_FEATURES, D_HIDDEN)
    b1 = np.random.randn(D_HIDDEN)
    W2 = np.random.randn(D_HIDDEN, D_HIDDEN)
    b2 = np.random.randn(D_HIDDEN)
    W3 = np.random.randn(D_HIDDEN, D_OUTPUT)
    b3 = np.random.randn(D_OUTPUT)

    # ReLU activation function
    def relu(x):
        return np.maximum(0, x)

    # --- 3. Forward Pass and Rank Calculation ---
    # Layer 1
    H1 = relu(X @ W1 + b1)
    rank_H1 = np.linalg.matrix_rank(H1)

    # Layer 2
    H2 = relu(H1 @ W2 + b2)
    rank_H2 = np.linalg.matrix_rank(H2)

    # Layer 3 (Output)
    H3 = relu(H2 @ W3 + b3)
    rank_H3 = np.linalg.matrix_rank(H3)

    # --- 4. Analysis and Output ---
    print("Analyzing the rank of latent representation matrices.")
    print("The rank of a matrix from a layer's output is limited by the rank of its input and its own dimensions.")
    print("\n--- Input ---")
    print(f"The rank of the input matrix X is {rank_X}.")

    print("\n--- Layer 1 Analysis ---")
    print(f"The rank of the output H1 is bounded by the rank of the input X.")
    print(f"Theoretical Bound: rank(H1) <= min(rank(X), dimension of layer 1)")
    print(f"rank(H1) <= min({rank_X}, {D_HIDDEN})")
    max_rank_H1 = min(rank_X, D_HIDDEN)
    print(f"So, the maximum possible rank for H1 is {max_rank_H1}.")
    print(f"Our simulation gives rank(H1) = {rank_H1}.")
    print("A. rank(H1) = 20: This is POSSIBLE because 20 is less than or equal to the max rank of 25.")
    print("B. rank(H1) = 50: This is IMPOSSIBLE because 50 is greater than the max rank of 25.")

    print("\n--- Layer 2 Analysis ---")
    print(f"The rank of H2 is bounded by the rank of its input, H1.")
    print(f"Theoretical Bound: rank(H2) <= min(rank(H1), dimension of layer 2)")
    print(f"rank(H2) <= min({max_rank_H1}, {D_HIDDEN})")
    max_rank_H2 = min(max_rank_H1, D_HIDDEN)
    print(f"So, the maximum possible rank for H2 is {max_rank_H2}.")
    print(f"Our simulation gives rank(H2) = {rank_H2}.")
    print("C. rank(H2) = 25: This is POSSIBLE, as it equals the maximum theoretical rank. This occurs if rank is preserved through the first two layers.")

    print("\n--- Layer 3 Analysis ---")
    print(f"The rank of H3 is bounded by its own dimensions and the rank of its input, H2.")
    print(f"Theoretical Bound: rank(H3) <= min(rank(H2), dimension of layer 3)")
    print(f"rank(H3) <= min({max_rank_H2}, {D_OUTPUT})")
    max_rank_H3 = min(max_rank_H2, D_OUTPUT)
    print(f"So, the maximum possible rank for H3 is {max_rank_H3}.")
    print(f"Our simulation gives rank(H3) = {rank_H3}.")
    print("D. rank(H3) = 15: This is IMPOSSIBLE because 15 is greater than the max rank of 10.")

    print("\n--- Final Conclusion ---")
    print("Statements that could be TRUE: A, C.")
    print("Statements that are FALSE: B, D.")
    print("The correct choice is H, which includes A and C.")

solve()
<<<H>>>