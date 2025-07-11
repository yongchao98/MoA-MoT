def analyze_network_ranks():
    """
    Analyzes and explains the possible ranks of latent representations
    in the described 3-layer MLP.
    """

    # --- Problem Definition ---
    n_datapoints = 100
    n_features = 50
    input_rank = 25

    layer_dims = {
        'layer1': 50,
        'layer2': 50,
        'layer3': 10
    }

    print("--- Rank Analysis Step-by-Step ---")
    print(f"Initial input matrix rank: {input_rank}")
    print("Key principles for rank calculation:")
    print("1. rank(A @ B) <= min(rank(A), rank(B))")
    print("2. rank(X @ W + b) <= rank(X @ W) + 1 (bias can increase rank by at most 1)")
    print("3. rank(ReLU(Z)) <= rank(Z) (ReLU is rank-non-increasing)")
    print("4. rank(M) <= min(rows, cols)")
    print("-" * 40)

    # --- Layer 1 Analysis ---
    print("\nAnalyzing Layer 1 (output H1):")
    # rank(X @ W1) <= min(rank(X), rank(W1))
    max_rank_after_linear = min(input_rank, layer_dims['layer1'])
    # rank(X@W1+b1) <= rank(X@W1) + 1
    max_rank_H1 = max_rank_after_linear + 1
    print(f"The maximum rank after the linear transform is min({input_rank}, {layer_dims['layer1']}) = {max_rank_after_linear}.")
    print(f"The maximum rank after adding bias and applying ReLU is {max_rank_after_linear} + 1 = {max_rank_H1}.")
    
    # Evaluate statements A and B
    stmt_A_rank = 20
    is_A_possible = stmt_A_rank <= max_rank_H1
    print(f"  > Statement A (rank={stmt_A_rank}): POSSIBLE, because {stmt_A_rank} <= {max_rank_H1}.")

    stmt_B_rank = 50
    is_B_possible = stmt_B_rank <= max_rank_H1
    print(f"  > Statement B (rank={stmt_B_rank}): IMPOSSIBLE, because {stmt_B_rank} > {max_rank_H1}.")
    print("-" * 40)

    # --- Layer 2 Analysis ---
    print("\nAnalyzing Layer 2 (output H2):")
    # The input to this layer is H1, with max rank max_rank_H1
    max_rank_after_linear_2 = min(max_rank_H1, layer_dims['layer2'])
    max_rank_H2 = max_rank_after_linear_2 + 1
    print(f"The input H1 has a maximum rank of {max_rank_H1}.")
    print(f"The maximum rank after the linear transform is min({max_rank_H1}, {layer_dims['layer2']}) = {max_rank_after_linear_2}.")
    print(f"The maximum rank after adding bias and applying ReLU is {max_rank_after_linear_2} + 1 = {max_rank_H2}.")

    # Evaluate statement C
    stmt_C_rank = 25
    is_C_possible = stmt_C_rank <= max_rank_H2
    print(f"  > Statement C (rank={stmt_C_rank}): POSSIBLE, because {stmt_C_rank} <= {max_rank_H2}.")
    print("-" * 40)

    # --- Layer 3 Analysis ---
    print("\nAnalyzing Layer 3 (output H3):")
    # The output dimension provides a strict upper bound on the rank.
    max_rank_H3 = layer_dims['layer3']
    print(f"The output matrix H3 has {layer_dims['layer3']} columns.")
    print(f"Therefore, its rank is strictly limited by its dimension, meaning the maximum possible rank is {max_rank_H3}.")

    # Evaluate statement D
    stmt_D_rank = 15
    is_D_possible = stmt_D_rank <= max_rank_H3
    print(f"  > Statement D (rank={stmt_D_rank}): IMPOSSIBLE, because {stmt_D_rank} > {max_rank_H3}.")
    print("-" * 40)

    # --- Conclusion ---
    print("\nConclusion:")
    print("The statements that could be true are A and C.")
    print("This corresponds to answer choice H.")

analyze_network_ranks()
print("\n<<<H>>>")