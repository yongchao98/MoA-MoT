import numpy as np
from scipy.stats import norm

def demonstrate_squared_norm_approximation(N):
    """
    Demonstrates that a shallow network with 2N hidden neurons can approximate
    the squared norm of an N-dimensional vector.

    Args:
        N (int): The dimension of the input vector.
    """
    if not isinstance(N, int) or N <= 0:
        print("N must be a positive integer.")
        return

    # 1. The GeLU activation function
    def gelu(x):
        # Gaussian Error Linear Unit
        return x * norm.cdf(x)

    # 2. Create a random N-dimensional input vector for testing.
    # We scale it down because our weight settings are derived from
    # Taylor series approximation around zero.
    x = np.random.randn(N) * 0.5
    
    # 3. The minimum hidden layer width is H = 2N
    H = 2 * N

    # --- Network Construction ---

    # 4. Hidden Layer (Input -> Hidden)
    # We set weights to be e_i and -e_i for i=1..N
    W1 = np.zeros((H, N))
    for i in range(N):
        W1[2 * i, i] = 1.0    # Corresponds to input x_i
        W1[2 * i + 1, i] = -1.0 # Corresponds to input -x_i
    
    # Hidden biases are all zero for this construction
    b1 = np.zeros(H)
    
    # Calculate hidden activations
    z1 = W1 @ x + b1
    a1 = gelu(z1)

    # 5. Output Layer (Hidden -> Output)
    # The output weights are set to a constant `c`. The value `sqrt(2*pi)` is
    # chosen to make the final output value numerically close to the target.
    c = np.sqrt(2 * np.pi)
    w2 = np.full(H, c)
    
    # Output bias is zero
    b2 = 0
    
    # Calculate pre-activation and final output
    z2 = w2 @ a1 + b2
    network_output = gelu(z2)

    # 6. Calculate the true squared norm
    true_value = np.sum(x**2)
    
    # --- Print Results ---
    print(f"Demonstration for an input vector in N={N} dimensions.")
    print(f"Input vector x: {np.round(x, 4)}")
    print("-" * 30)
    print(f"Target (True squared norm): {true_value:.6f}")
    print(f"Network Prediction (H={H}): {network_output:.6f}")
    print(f"Approximation Error: {abs(true_value - network_output):.6f}")
    print("-" * 30)

    # Output the final equation for the minimum width
    print("The minimum hidden-layer width (H) is given by the formula in terms of N:")
    print(f"H = 2 * {N} = {2 * N}")

# Run the demonstration for a chosen N
# You can change this value to test other dimensions.
example_N = 3
demonstrate_squared_norm_approximation(example_N)