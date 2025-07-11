import numpy as np

def analyze_network_rank():
    """
    Analyzes and demonstrates the rank of latent representations in a 3-layer MLP.
    """
    # Define network and data parameters
    n_samples = 100
    n_features = 50
    input_rank = 25
    layer1_dim = 50
    layer2_dim = 50
    layer3_dim = 10

    # For reproducibility
    np.random.seed(42)

    # 1. Create an input matrix X with shape (100, 50) and rank 25
    # Create two random matrices U (100x25) and V (25x50)
    U = np.random.randn(n_samples, input_rank)
    V = np.random.randn(input_rank, n_features)
    # Their product X will have a rank of at most 25
    X = U @ V
    rank_X = np.linalg.matrix_rank(X)

    # 2. Define the 3-layer MLP
    # Layer 1 weights and biases
    W1 = np.random.randn(n_features, layer1_dim)
    b1 = np.random.randn(1, layer1_dim)
    # Layer 2 weights and biases
    W2 = np.random.randn(layer1_dim, layer2_dim)
    b2 = np.random.randn(1, layer2_dim)
    # Layer 3 weights and biases
    W3 = np.random.randn(layer2_dim, layer3_dim)
    b3 = np.random.randn(1, layer3_dim)

    # ReLU activation function
    relu = lambda x: np.maximum(0, x)

    # 3. Forward pass and calculate ranks
    # Layer 1
    A1 = relu(X @ W1 + b1)
    rank_A1 = np.linalg.matrix_rank(A1)
    # Layer 2
    A2 = relu(A1 @ W2 + b2)
    rank_A2 = np.linalg.matrix_rank(A2)
    # Layer 3
    A3 = relu(A2 @ W3 + b3)
    rank_A3 = np.linalg.matrix_rank(A3)

    # 4. Output the theoretical analysis and simulation results
    print("--- Theoretical Rank Analysis ---")
    print(f"Input matrix X has shape ({n_samples}, {n_features}) and rank = {input_rank}")

    # Analysis for Layer 1
    max_rank_A1 = min(rank_X + 1, layer1_dim)
    print("\nLayer 1 (output shape 100x50):")
    print(f"The rank of the latent representation is bounded by min(rank(X) + 1, dim) = min({rank_X} + 1, {layer1_dim}) = {max_rank_A1}.")
    print(f"  A. rank=20: Possible, because 20 <= {max_rank_A1}.")
    print(f"  B. rank=50: Impossible, because 50 > {max_rank_A1}.")

    # Analysis for Layer 2
    max_rank_A2 = min(max_rank_A1 + 1, layer2_dim)
    print("\nLayer 2 (output shape 100x50):")
    print(f"The rank is bounded by the rank of its input (max {max_rank_A1}). So its rank is at most min({max_rank_A1} + 1, {layer2_dim}) = {max_rank_A2}.")
    print(f"  C. rank=25: Possible, because 25 <= {max_rank_A2}.")

    # Analysis for Layer 3
    max_rank_A3 = min(n_samples, layer3_dim)
    print("\nLayer 3 (output shape 100x10):")
    print(f"The rank is bounded by its dimension: min(100, 10) = {max_rank_A3}.")
    print(f"  D. rank=15: Impossible, because 15 > {max_rank_A3}.")

    print("\n--- Rank Simulation (with one random initialization) ---")
    print(f"Actual rank of the generated input matrix X: {rank_X}")
    print(f"Actual rank of Layer 1 latent representation (A1): {rank_A1}")
    print(f"Actual rank of Layer 2 latent representation (A2): {rank_A2}")
    print(f"Actual rank of Layer 3 latent representation (A3): {rank_A3}")
    
    print("\nConclusion: Statements A and C are the only ones that could be true.")

if __name__ == '__main__':
    analyze_network_rank()
<<<H>>>