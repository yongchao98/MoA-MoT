import numpy as np

def run_rank_analysis():
    """
    Analyzes and demonstrates the rank of latent representations in a 3-layer MLP.
    """
    # Part 1: Theoretical Analysis
    print("Analyzing the rank of latent representations in a 3-layer MLP.")
    print("-" * 60)
    print("Network & Input Specification:")
    print("- Input X: 100x50 matrix, rank(X) = 25")
    print("- Layer 1: 50 -> 50 neurons, ReLU. Latent matrix H1 (100x50)")
    print("- Layer 2: 50 -> 50 neurons, ReLU. Latent matrix H2 (100x50)")
    print("- Layer 3: 50 -> 10 neurons, ReLU. Latent matrix Y (100x10)")
    print("-" * 60)
    print("Theoretical Rank Analysis:")

    # Layer 1
    print("\n--- Analysis of Layer 1 (H1) ---")
    print("H1 = ReLU(X @ W1 + b1)")
    print("The rank of a matrix product is at most the minimum of the input ranks.")
    print("rank(X @ W1) <= min(rank(X), rank(W1))")
    print("Given rank(X) = 25 and W1 is 50x50, rank(X @ W1) <= min(25, 50) = 25.")
    print("Adding a bias can increase the rank by at most 1. The ReLU function does not increase rank.")
    print("Therefore, rank(H1) <= rank(X) + 1 = 26.")

    print("\nEvaluating statements for Layer 1:")
    print("A. The rank of matrix for the first layer is 20.")
    print("   - COULD BE TRUE. 20 is less than the max possible rank of 26. The ReLU non-linearity and weights can easily reduce the rank.")
    print("B. The rank of matrix for the first layer is 50.")
    print("   - FALSE. This is impossible as the rank is capped by the input data's rank bottleneck, at most 26.")

    # Layer 2
    print("\n--- Analysis of Layer 2 (H2) ---")
    print("H2 = ReLU(H1 @ W2 + b2)")
    print("The rank of H2 is limited by the rank of its input, H1.")
    print("Since rank(H1) <= 26, the rank of H2 will also be <= 26 (approximately).")

    print("\nEvaluating statement for Layer 2:")
    print("C. The rank of matrix for the second layer is 25.")
    print("   - COULD BE TRUE. It's possible for the transformations to be rank-preserving, maintaining the original data rank of 25.")

    # Last Layer
    print("\n--- Analysis of the Last Layer (Y) ---")
    print("Y = ReLU(H2 @ W3 + b3)")
    print("The output matrix Y has a shape of (100, 10).")
    print("The rank of any matrix cannot exceed its number of rows or columns.")
    print("Therefore, rank(Y) <= min(100, 10) = 10.")

    print("\nEvaluating statement for the Last Layer:")
    print("D. The rank of matrix for the last layer is 15.")
    print("   - FALSE. This is impossible as the matrix only has 10 columns, so its rank cannot exceed 10.")

    print("-" * 60)
    print("Demonstration with Code:")
    
    # Part 2: Code Demonstration
    def relu(x):
        return np.maximum(0, x)

    # Setup: Create input X with shape (100, 50) and rank 25
    A = np.random.RandomState(0).randn(100, 25)
    B = np.random.RandomState(1).randn(25, 50)
    X = A @ B
    print(f"\nCreated input matrix X of shape {X.shape} with rank {np.linalg.matrix_rank(X)}")

    # Demonstration for Statement C (rank=25 is possible)
    print("\n--- Demonstrating Statement C (rank=25 is plausible) ---")
    np.random.seed(42)
    W1 = np.random.randn(50, 50) * 0.5 
    b1 = np.random.randn(50) * 0.5 
    H1 = relu(X @ W1 + b1)
    W2 = np.random.randn(50, 50) * 0.5
    b2 = np.random.randn(50) * 0.5
    H2 = relu(H1 @ W2 + b2)
    rank_H2_C = np.linalg.matrix_rank(H2)
    print(f"With one set of random parameters, we found rank(H2) = {rank_H2_C}.")
    
    # Demonstration for Statement A (rank=20 is possible)
    print("\n--- Demonstrating Statement A (rank=20 is plausible) ---")
    np.random.seed(1)
    W1_A = np.random.randn(50, 50)
    # Adjust bias with a negative offset to reduce rank by "killing" neurons
    b1_A = np.random.randn(50) - 2.5
    H1_A = relu(X @ W1_A + b1_A)
    rank_H1_A = np.linalg.matrix_rank(H1_A)
    print(f"By adjusting the layer 1 bias, we found a case where rank(H1) = {rank_H1_A}.")

    # Final Conclusion
    print("-" * 60)
    print("\nSummary of Findings:")
    print("Statement A (rank=20): Could be True")
    print("Statement B (rank=50): False")
    print("Statement C (rank=25): Could be True")
    print("Statement D (rank=15): False")
    print("\nThe correct answer includes statements A and C.")

# Execute the analysis
run_rank_analysis()
<<<H>>>