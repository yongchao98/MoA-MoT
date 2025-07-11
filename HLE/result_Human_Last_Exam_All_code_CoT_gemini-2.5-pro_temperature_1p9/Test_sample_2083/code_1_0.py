import numpy as np

# This script demonstrates that a shallow neural network with a hidden-layer
# width of 2*N can approximate the squared norm of an N-dimensional vector.
# The explanation for why 2*N is the *minimum* width is provided in the text.

def gelu_approx(x):
    """
    An accurate approximation of the GeLU function.
    GeLU(x) = 0.5 * x * (1 + tanh(sqrt(2/pi) * (x + 0.044715 * x^3)))
    """
    return 0.5 * x * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x**3)))

def demonstrate_squared_norm_computation(N, input_vector):
    """
    Constructs a 1-hidden-layer network to compute ||x||^2 and shows it works.
    N: int, dimensionality of the input vector
    input_vector: np.array, the input vector of shape (N,)
    """
    print(f"--- Demonstration for N = {N} ---")
    print(f"Input vector x = {input_vector}")

    # --- Network Construction ---
    # The minimum required width is 2*N
    H = 2 * N
    print(f"Using a hidden layer width of H = 2 * N = {H}")

    # 1. Input to Hidden Layer Weights and Biases
    # W1 is an H x N matrix, b1 is a vector of size H
    W1 = np.zeros((H, N))
    for i in range(N):
        W1[2 * i, i] = 1.0      # Neuron 2*i is fed x_i
        W1[2 * i + 1, i] = -1.0 # Neuron 2*i+1 is fed -x_i
    
    b1 = np.zeros(H) # Biases are all zero

    # 2. Hidden Layer to Output Weights and Biases
    # w2 is a vector of size H, b2 is a scalar
    # From our derivation, GeLU(z) + GeLU(-z) ≈ sqrt(2/pi) * z^2
    # So, z^2 ≈ sqrt(pi/2) * (GeLU(z) + GeLU(-z))
    # The output weight C for each neuron must be sqrt(pi/2).
    C = np.sqrt(np.pi / 2.0)
    w2 = np.full(H, C)
    b2 = 0.0 # Output bias is zero

    # --- Forward Pass ---
    # hidden_activations = gelu(W1 @ x + b1)
    hidden_layer_input = W1 @ input_vector + b1
    hidden_layer_activations = gelu_approx(hidden_layer_input)
    
    # network_output = w2 @ hidden_activations + b2
    network_output = w2 @ hidden_layer_activations + b2

    # --- Compare with the true value ---
    true_value = np.sum(input_vector**2)

    print(f"Network output (approximation of ||x||^2): {network_output:.6f}")
    print(f"True value of ||x||^2: {true_value:.6f}")
    
    # --- Final Conclusion ---
    # The final equation gives the minimum width in terms of N.
    print("\nBased on the analysis, the minimum required width is given by the equation:")
    # We output each component of the final equation as requested.
    number_of_dimensions = "N"
    coefficient = 2
    print(f"Width = {coefficient} * {number_of_dimensions}")


# --- Main Execution ---
# Set the dimension N for the demonstration
N = 4

# Create a random input vector. Using smaller values demonstrates the
# local approximation principle more clearly.
np.random.seed(42)
x = np.random.randn(N) * 0.5 

demonstrate_squared_norm_computation(N, x)