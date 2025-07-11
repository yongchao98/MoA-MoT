import numpy as np
from scipy.special import erf

def solve():
    """
    This script determines the minimum hidden-layer width for a shallow GeLU
    network to compute the squared norm of an N-dimensional vector and
    demonstrates the result.
    """
    # N is the dimension of the input vector. You can change this value.
    N = 3

    # According to approximation theory for neural networks, the minimum
    # hidden-layer width H required is 2*N.
    H = 2 * N

    # Print the final answer equation
    print(f"For an input vector of dimension N = {N}")
    print(f"The minimum required hidden layer width (H) is given by the equation:")
    print(f"H = 2 * N")
    print(f"H = 2 * {N} = {H}")
    print("-" * 40)
    print("Demonstration of the principle:")

    # 1. Define the GeLU activation function
    def gelu(x):
        return 0.5 * x * (1 + erf(x / np.sqrt(2)))

    # 2. Define network parameters based on the 2N construction.
    # The pre-output A(x) is constructed to approximate 2 * ||x||^2.
    # The Taylor expansion of gelu(z)+gelu(-z) is sqrt(2/pi)*z^2 + O(z^4).
    # To make A(x) = v * sum(gelu(x_i)+gelu(-x_i)) approx 2*||x||^2,
    # we need v * sqrt(2/pi) = 2, which means v = sqrt(2*pi).
    output_weight_val = np.sqrt(2 * np.pi)

    # Hidden layer weights (H x N matrix)
    # For each dimension i, we use weights e_i and -e_i.
    hidden_weights = np.zeros((H, N))
    for i in range(N):
        hidden_weights[2*i, i] = 1.0
        hidden_weights[2*i + 1, i] = -1.0

    # Hidden layer biases are all zero for this construction
    hidden_biases = np.zeros(H)

    # Output layer weights are all equal
    output_weights = np.full(H, output_weight_val)

    # Output bias is zero
    output_bias = 0.0

    # 3. Generate a random input vector x.
    # The approximation is most accurate for small x where GeLU(z) ~ z/2.
    np.random.seed(42)
    x = np.random.randn(N) * 0.5

    # 4. Compute the network's output
    # Hidden layer pre-activation: z_k = w_k^T * x + b_k
    z_hidden = hidden_weights @ x + hidden_biases
    # Hidden layer activation: h_k = GeLU(z_k)
    h_hidden = gelu(z_hidden)
    # Output layer pre-activation: A(x) = v^T * h + c
    z_output = np.dot(output_weights, h_hidden) + output_bias
    # Final network output: y = GeLU(A(x))
    y_network = gelu(z_output)

    # 5. Compute the true value for comparison
    y_true = np.sum(x**2)

    # 6. Print the results
    print(f"Random input vector x:          {np.round(x, 4)}")
    print(f"True squared norm ||x||^2:      {y_true:.6f}")
    print(f"Network (H={H}) output y(x):   {y_network:.6f}")
    print(f"Approximation Error:          {abs(y_true - y_network):.6f}")

solve()