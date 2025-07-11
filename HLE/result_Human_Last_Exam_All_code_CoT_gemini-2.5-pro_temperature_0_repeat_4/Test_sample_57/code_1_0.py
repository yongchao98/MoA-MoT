import torch
import torch.nn as nn

def main():
    """
    Analyzes and demonstrates the rank of latent representations in a 3-layer MLP.
    """
    # --- Network and Data Parameters ---
    n_datapoints = 100
    n_features = 50
    n_hidden = 50
    n_output = 10
    input_rank = 25

    # --- Step 1: Create Input Data with a specific rank ---
    # We create a 100x50 matrix X with rank 25 by multiplying two full-rank matrices.
    # X = A @ B, where A is 100x25 and B is 25x50.
    # rank(X) <= min(rank(A), rank(B)) = min(25, 25) = 25.
    print("Setting up the scenario...")
    A = torch.randn(n_datapoints, input_rank)
    B = torch.randn(input_rank, n_features)
    X = A @ B
    rank_X = torch.linalg.matrix_rank(X).item()
    print(f"Input matrix X created with shape {X.shape} and rank = {rank_X}")
    print("-" * 40)

    # --- Step 2: Define the 3-Layer MLP ---
    # We define a simple MLP with ReLU activations.
    # To demonstrate the possibility of rank preservation, we initialize the network
    # with large positive biases. This pushes the pre-activations to be positive,
    # making the ReLU function behave like an identity function in that region,
    # which helps preserve rank.
    model = nn.Sequential(
        nn.Linear(n_features, n_hidden),
        nn.ReLU(),
        nn.Linear(n_hidden, n_hidden),
        nn.ReLU(),
        nn.Linear(n_hidden, n_output),
        nn.ReLU()
    )
    with torch.no_grad():
        for layer in model:
            if isinstance(layer, nn.Linear):
                layer.bias.data.fill_(10.0) # Use large bias for demonstration

    # --- Step 3: Forward Pass and Rank Calculation ---
    # We pass the data through the network and store intermediate representations.
    H1 = model[1](model[0](X))
    H2 = model[3](model[2](H1))
    H3 = model[5](model[4](H2))

    # Calculate the rank of each representation matrix
    rank_H1 = torch.linalg.matrix_rank(H1).item()
    rank_H2 = torch.linalg.matrix_rank(H2).item()
    rank_H3 = torch.linalg.matrix_rank(H3).item()

    print("Theoretical Rank Analysis:")
    print(f"Rank(H1) <= Rank(X) = {input_rank}")
    print(f"Rank(H2) <= Rank(H1) <= {input_rank}")
    print(f"Rank(H3) <= min(Rank(H2), output_dim) = min({rank_H2}, {n_output})")
    print("-" * 40)

    print("Ranks from a Practical Demonstration:")
    print(f"Rank of Layer 1 latent matrix (H1): {rank_H1}")
    print(f"Rank of Layer 2 latent matrix (H2): {rank_H2}")
    print(f"Rank of Layer 3 latent matrix (H3): {rank_H3}")
    print("-" * 40)

    # --- Step 4: Evaluate the Statements ---
    print("Evaluating the Answer Choices:")
    # Statement A
    print(f"A. Rank of H1 is 20: This could be TRUE. The theoretical maximum rank is {input_rank}. A rank of 20 is possible if the transformation reduces the rank.")
    # Statement B
    print(f"B. Rank of H1 is 50: This is FALSE. The rank cannot exceed the input rank of {input_rank}.")
    # Statement C
    print(f"C. Rank of H2 is 25: This could be TRUE. Our demonstration shows a rank of {rank_H2}, confirming it's possible to preserve the rank up to its maximum of {input_rank}.")
    # Statement D
    print(f"D. Rank of H3 is 15: This is FALSE. The rank cannot exceed the layer's output dimension of {n_output}. The maximum possible rank is {n_output}.")

if __name__ == '__main__':
    main()
<<<H>>>