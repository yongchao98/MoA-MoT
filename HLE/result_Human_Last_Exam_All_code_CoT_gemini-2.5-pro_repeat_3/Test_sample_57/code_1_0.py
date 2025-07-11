import torch

def solve():
    """
    Analyzes and demonstrates the rank of latent representations in an MLP.
    """
    # Network and data parameters
    n_samples = 100
    n_features = 50
    n_hidden = 50
    n_output = 10
    rank_input = 25

    # 1. Create an input matrix X with a specified rank
    # To create a matrix of a specific rank, we can multiply two smaller matrices.
    # rank(A @ B) <= min(rank(A), rank(B))
    A = torch.randn(n_samples, rank_input)
    B = torch.randn(rank_input, n_features)
    X = A @ B
    # For random A and B, rank(X) will be rank_input with very high probability.
    rank_X = torch.linalg.matrix_rank(X).item()

    # 2. Define the 3-layer MLP layers (weights and biases)
    # With random initialization, weight matrices will be full-rank.
    W1 = torch.randn(n_features, n_hidden)
    b1 = torch.randn(n_hidden)
    W2 = torch.randn(n_hidden, n_hidden)
    b2 = torch.randn(n_hidden)
    W3 = torch.randn(n_hidden, n_output)
    b3 = torch.randn(n_output)

    # 3. Perform the forward pass and compute latent representations
    # Layer 1
    H1 = torch.relu(X @ W1 + b1)
    # Layer 2
    H2 = torch.relu(H1 @ W2 + b2)
    # Layer 3
    H3 = torch.relu(H2 @ W3 + b3)

    # 4. Calculate the rank of each latent representation matrix
    rank_H1 = torch.linalg.matrix_rank(H1).item()
    rank_H2 = torch.linalg.matrix_rank(H2).item()
    rank_H3 = torch.linalg.matrix_rank(H3).item()

    print("--- Rank Analysis ---")
    print(f"The rank of the input matrix X is: {rank_X}")
    print(f"The rank of the first layer's latent representation H1 is: {rank_H1}")
    print(f"The rank of the second layer's latent representation H2 is: {rank_H2}")
    print(f"The rank of the last layer's latent representation H3 is: {rank_H3}")
    print("\n--- Theoretical Evaluation of Statements ---")

    # Statement A
    print("Statement A: The rank of matrix containing latent space representations of the first layer is 20.")
    print(f"Theoretical Bound: rank(H1) <= rank(X) = {rank_input}.")
    print("Conclusion: This is POSSIBLE. A rank of 20 is less than or equal to the theoretical maximum of 25.\n")

    # Statement B
    print("Statement B: The rank of matrix containing latent space representations of the first layer is 50.")
    print(f"Theoretical Bound: rank(H1) <= rank(X) = {rank_input}.")
    print(f"Conclusion: This is IMPOSSIBLE. A rank of 50 exceeds the theoretical maximum of 25.\n")

    # Statement C
    print("Statement C: The rank of matrix containing latent space representations of the second layer is 25.")
    print(f"Theoretical Bound: rank(H2) <= rank(H1) <= {rank_input}.")
    print("Conclusion: This is POSSIBLE. It is plausible for the rank to be maintained at its maximum through the first two layers.\n")

    # Statement D
    print("Statement D: The rank of matrix containing latent space representations of the last layer is 15.")
    print(f"Theoretical Bound: rank(H3) <= min(rows, columns) = min({n_samples}, {n_output}) = {n_output}.")
    print(f"Conclusion: This is IMPOSSIBLE. A rank of 15 exceeds the layer's output dimension of 10.\n")

solve()
<<<H>>>