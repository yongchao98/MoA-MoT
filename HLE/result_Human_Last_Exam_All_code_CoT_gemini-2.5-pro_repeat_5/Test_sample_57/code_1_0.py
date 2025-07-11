def analyze_network_rank():
    """
    Analyzes and explains the rank of latent representations in a 3-layer MLP.
    """
    # Network and Input descriptions
    n_samples = 100
    n_features = 50
    initial_rank = 25

    layer1_dim = 50
    layer2_dim = 50
    layer3_dim = 10

    print("Analyzing the rank of latent representation matrices.\n")
    print(f"Input matrix X has shape ({n_samples}, {n_features}) and rank = {initial_rank}.\n")

    # --- Layer 1 Analysis ---
    print("--- Layer 1 Analysis ---")
    # H1 = ReLU(X @ W1 + b1) where W1 is 50x50
    # The rank of W1 can be at most 50. We assume it is full rank for max rank calculation.
    rank_W1 = 50
    # Step 1: Linear transformation (X @ W1)
    # rank(X @ W1) <= min(rank(X), rank(W1))
    max_rank_after_linear1 = min(initial_rank, rank_W1)
    print(f"The rank after the first linear transformation (X @ W1) is at most:")
    print(f"rank(X @ W1) <= min(rank(X), rank(W1)) = min({initial_rank}, {rank_W1}) = {max_rank_after_linear1}")

    # Step 2: Bias addition (+ b1)
    # rank(X @ W1 + b1) <= rank(X @ W1) + 1
    max_rank_after_bias1 = max_rank_after_linear1 + 1
    print(f"\nAdding the bias term can increase the rank by at most 1:")
    print(f"rank(X @ W1 + b1) <= rank(X @ W1) + 1 <= {max_rank_after_linear1} + 1 = {max_rank_after_bias1}")

    # Step 3: ReLU activation
    # rank(H1) = rank(ReLU(...)) <= rank(...)
    max_rank_H1 = max_rank_after_bias1
    print(f"\nThe ReLU activation cannot increase the rank, so the final maximum rank for the first layer's latent representation H1 is:")
    print(f"rank(H1) <= {max_rank_H1}")

    print("\nEvaluating statements for Layer 1:")
    # A. The rank of matrix containing latent space representations of the first layer is 20.
    rank_A = 20
    print(f"  A. rank(H1) = {rank_A}: This is possible because {rank_A} <= {max_rank_H1}.")
    # B. The rank of matrix containing latent space representations of the first layer is 50.
    rank_B = 50
    print(f"  B. rank(H1) = {rank_B}: This is FALSE because {rank_B} > {max_rank_H1}.\n")


    # --- Layer 2 Analysis ---
    print("--- Layer 2 Analysis ---")
    # H2 = ReLU(H1 @ W2 + b2) where W2 is 50x50
    # The input to this layer is H1, with a maximum possible rank of max_rank_H1.
    rank_W2 = 50
    # Step 1: Linear transformation (H1 @ W2)
    max_rank_after_linear2 = min(max_rank_H1, rank_W2)
    print(f"The rank after the second linear transformation (H1 @ W2) is at most:")
    print(f"rank(H1 @ W2) <= min(rank(H1), rank(W2)) <= min({max_rank_H1}, {rank_W2}) = {max_rank_after_linear2}")

    # Step 2: Bias addition (+ b2)
    max_rank_after_bias2 = max_rank_after_linear2 + 1
    print(f"\nAdding the bias term can increase the rank by at most 1:")
    print(f"rank(H1 @ W2 + b2) <= rank(H1 @ W2) + 1 <= {max_rank_after_linear2} + 1 = {max_rank_after_bias2}")

    # Step 3: ReLU activation
    max_rank_H2 = max_rank_after_bias2
    print(f"\nThe ReLU activation cannot increase the rank, so the final maximum rank for the second layer's latent representation H2 is:")
    print(f"rank(H2) <= {max_rank_H2}")

    print("\nEvaluating statement for Layer 2:")
    # C. The rank of matrix containing latent space representations of the second layer is 25.
    rank_C = 25
    print(f"  C. rank(H2) = {rank_C}: This is possible because {rank_C} <= {max_rank_H2}.\n")


    # --- Layer 3 Analysis ---
    print("--- Layer 3 Analysis ---")
    # H3 = ReLU(H2 @ W3 + b3) where H3 is 100x10
    # The rank is fundamentally limited by the matrix dimensions.
    max_rank_H3_by_dim = min(n_samples, layer3_dim)
    print(f"The output of the last layer, H3, has a shape of ({n_samples}, {layer3_dim}).")
    print("The rank of a matrix cannot be larger than its number of rows or columns.")
    print(f"Therefore, the maximum possible rank for H3 is constrained by its dimensions:")
    print(f"rank(H3) <= min({n_samples}, {layer3_dim}) = {max_rank_H3_by_dim}")

    print("\nEvaluating statement for Layer 3:")
    # D. The rank of matrix containing latent space representations of the last layer is 15.
    rank_D = 15
    print(f"  D. rank(H3) = {rank_D}: This is FALSE because {rank_D} > {max_rank_H3_by_dim}.\n")


    # --- Conclusion ---
    print("--- Conclusion ---")
    print("The statements that could be true are A and C.")
    print("A: rank(H1)=20 (Possible, since 20 <= 26)")
    print("B: rank(H1)=50 (Impossible, since 50 > 26)")
    print("C: rank(H2)=25 (Possible, since 25 <= 27)")
    print("D: rank(H3)=15 (Impossible, since 15 > 10)")
    print("\nTherefore, the correct choice is H, which includes A and C.")

if __name__ == '__main__':
    analyze_network_rank()