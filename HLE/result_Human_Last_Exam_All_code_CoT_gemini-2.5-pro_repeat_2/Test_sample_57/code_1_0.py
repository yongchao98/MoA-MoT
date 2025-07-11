import torch

def solve_and_explain():
    """
    Analyzes and demonstrates the rank of latent representations in a MLP.
    """
    # --- 1. Theoretical Analysis based on Rank Properties ---
    print("--- Theoretical Analysis of Matrix Rank Propagation ---")

    # Input Layer
    rank_X = 25
    print(f"The problem states the input matrix X has a rank of {rank_X}.")
    print("A linear transformation (like multiplying by weights) or a non-linear activation (like ReLU) cannot increase the rank of a matrix.")
    print("-" * 50)

    # --- Analysis for Layer 1 ---
    layer1_output_dim = 50
    # The rank of the first latent representation, L1, is limited by the rank of the input, X.
    # Equation: rank(L1) <= rank(X)
    max_rank_L1 = rank_X
    print(f"Layer 1 (L1) Analysis:")
    print(f"The rank of L1 is limited by the rank of the input X. So, rank(L1) <= {max_rank_L1}.")
    # Evaluate statements A and B
    print("  - Statement A: 'The rank of matrix ... of the first layer is 20.'")
    print(f"    This COULD BE TRUE, because 20 is less than or equal to the maximum possible rank of {max_rank_L1}.")
    print("  - Statement B: 'The rank of matrix ... of the first layer is 50.'")
    print(f"    This IS IMPOSSIBLE, because 50 is greater than the maximum possible rank of {max_rank_L1}.")
    print("-" * 50)

    # --- Analysis for Layer 2 ---
    layer2_output_dim = 50
    # The rank of L2 is limited by the rank of its input, L1.
    # Equation: rank(L2) <= rank(L1) <= rank(X)
    max_rank_L2 = rank_X
    print(f"Layer 2 (L2) Analysis:")
    print(f"The rank of L2 is limited by the rank of L1. Since rank(L1) <= {max_rank_L1}, it follows that rank(L2) <= {max_rank_L2}.")
    # Evaluate statement C
    print("  - Statement C: 'The rank of matrix ... of the second layer is 25.'")
    print(f"    This COULD BE TRUE. It is possible for the rank to be maintained at {max_rank_L2} through the first two layers.")
    print("-" * 50)

    # --- Analysis for the Last Layer ---
    last_layer_output_dim = 10
    # The rank of L3 is limited by its own dimensions, specifically the smaller one (10).
    # Equation: rank(L3) <= min(number of rows, number of columns)
    max_rank_L3 = last_layer_output_dim
    print(f"Last Layer (L3) Analysis:")
    print(f"The rank of L3 is limited by its output dimension, which is {last_layer_output_dim}.")
    print(f"Regardless of the input rank, rank(L3) must be <= {max_rank_L3}.")
    # Evaluate statement D
    print("  - Statement D: 'The rank of matrix ... of the last layer is 15.'")
    print(f"    This IS IMPOSSIBLE, because 15 is greater than the maximum possible rank of {max_rank_L3}.")
    print("-" * 50)

    # --- Conclusion ---
    print("\nConclusion:")
    print("The statements that could be true are A and C.")
    print("\n--- Numerical Simulation to Verify ---")
    print("The following is a simulation showing one possible outcome consistent with our theory.")
    
    # Set seed for reproducibility
    torch.manual_seed(42)

    # Create an input matrix X with the specified rank
    A = torch.randn(100, rank_X)
    B = torch.randn(rank_X, 50)
    X = A @ B

    # Define network layers
    W1 = torch.randn(50, layer1_output_dim)
    b1 = torch.randn(layer1_output_dim)
    W2 = torch.randn(layer1_output_dim, layer2_output_dim)
    b2 = torch.randn(layer2_output_dim)
    W3 = torch.randn(layer2_output_dim, last_layer_output_dim)
    b3 = torch.randn(last_layer_output_dim)
    relu = torch.nn.ReLU()

    # Forward pass
    L1 = relu(X @ W1 + b1)
    L2 = relu(L1 @ W2 + b2)
    L3 = relu(L2 @ W3 + b3)

    # Calculate and print the ranks from the simulation
    print(f"Simulated rank of X:  {torch.linalg.matrix_rank(X).item()}")
    print(f"Simulated rank of L1: {torch.linalg.matrix_rank(L1).item()} (Theoretical max: {max_rank_L1})")
    print(f"Simulated rank of L2: {torch.linalg.matrix_rank(L2).item()} (Theoretical max: {max_rank_L2})")
    print(f"Simulated rank of L3: {torch.linalg.matrix_rank(L3).item()} (Theoretical max: {max_rank_L3})")

solve_and_explain()
<<<H>>>