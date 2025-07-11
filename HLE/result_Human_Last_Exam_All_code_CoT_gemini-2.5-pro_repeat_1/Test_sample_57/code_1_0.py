import numpy as np

def analyze_network_rank():
    """
    Analyzes the maximum possible rank of latent representations in a 3-layer MLP
    based on the problem description.
    """
    # Network and Input Description
    n_samples = 100
    n_features = 50
    input_rank = 25

    dim_layer1 = 50
    dim_layer2 = 50
    dim_output = 10

    print("Analyzing the rank of latent representation matrices...\n")

    # --- Layer 1 Analysis ---
    # The rank after a linear transformation is at most the rank of the input.
    # rank(X @ W1) <= min(rank(X), rank(W1))
    # Adding a bias term can increase the rank by at most 1.
    # The ReLU activation can preserve or decrease the rank, but not increase it.
    max_rank_h1 = min(input_rank + 1, dim_layer1)
    
    print("--- Layer 1 ---")
    print(f"The input matrix X has rank {input_rank}.")
    print(f"After the linear transformation and bias addition, the rank is at most {input_rank} + 1 = {input_rank + 1}.")
    print(f"The ReLU activation cannot increase the rank, so the maximum possible rank for the first latent representation (H1) is {max_rank_h1}.")
    
    # Evaluate statements A and B
    rank_A = 20
    is_A_possible = rank_A <= max_rank_h1
    print(f"\nStatement A: The rank of H1 is {rank_A}.")
    print(f"Analysis: {rank_A} <= {max_rank_h1} is {is_A_possible}. This statement could be True.")

    rank_B = 50
    is_B_possible = rank_B <= max_rank_h1
    print(f"\nStatement B: The rank of H1 is {rank_B}.")
    print(f"Analysis: {rank_B} <= {max_rank_h1} is {is_B_possible}. This statement is False.")

    # --- Layer 2 Analysis ---
    # The input to this layer is H1, with a maximum rank of max_rank_h1.
    max_rank_h2 = min(max_rank_h1 + 1, dim_layer2)
    
    print("\n--- Layer 2 ---")
    print(f"The input to this layer (H1) has a maximum rank of {max_rank_h1}.")
    print(f"Following the same logic, the maximum possible rank for the second latent representation (H2) is {max_rank_h2}.")
    
    # Evaluate statement C
    rank_C = 25
    is_C_possible = rank_C <= max_rank_h2
    print(f"\nStatement C: The rank of H2 is {rank_C}.")
    print(f"Analysis: {rank_C} <= {max_rank_h2} is {is_C_possible}. This statement could be True.")

    # --- Last Layer Analysis ---
    # The output dimension provides the tightest constraint here.
    # The output matrix H3 has shape (100, 10), so its rank cannot exceed 10.
    max_rank_h3 = dim_output

    print("\n--- Last Layer ---")
    print(f"The output of the last layer (H3) is a matrix of size ({n_samples}, {dim_output}).")
    print(f"The rank of a matrix cannot exceed its smallest dimension.")
    print(f"Therefore, the maximum possible rank for the last latent representation (H3) is {max_rank_h3}.")

    # Evaluate statement D
    rank_D = 15
    is_D_possible = rank_D <= max_rank_h3
    print(f"\nStatement D: The rank of H3 is {rank_D}.")
    print(f"Analysis: {rank_D} <= {max_rank_h3} is {is_D_possible}. This statement is False.")
    
    print("\n--- Conclusion ---")
    print("Statements that could be True are A and C.")

if __name__ == '__main__':
    analyze_network_rank()
    print("\n<<<H>>>")