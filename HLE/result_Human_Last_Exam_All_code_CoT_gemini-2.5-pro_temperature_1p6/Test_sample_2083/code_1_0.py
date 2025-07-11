import numpy as np
from scipy.stats import norm

def gelu(x):
    """Gaussian Error Linear Unit (GELU) activation function."""
    return x * norm.cdf(x)

def squared_norm_network(x, alpha=0.1):
    """
    Computes the squared norm of a vector x using a shallow neural network
    with 2N hidden neurons and GELU activation.

    Args:
        x (np.ndarray): The N-dimensional input vector.
        alpha (float): A scaling parameter. Smaller values give better
                       accuracy for a larger range of x but might lead to
                       vanishing outputs if not compensated.

    Returns:
        float: The approximated squared norm of x.
    """
    N = len(x)
    H = 2 * N

    # The Taylor series for g(z)+g(-z) around z=0 is (2/sqrt(2*pi)) * z^2 + O(z^4)
    # Let z = alpha * x_j. We want the coefficient of x_j^2 to be 1.
    # The output of the j-th pair is approx C * (2/sqrt(2*pi)) * (alpha*x_j)^2
    # We set C * (2/sqrt(2*pi)) * alpha^2 = 1 => C = sqrt(2*pi) / (2 * alpha^2)
    output_weight = np.sqrt(2 * np.pi) / (2 * alpha**2)

    hidden_layer_sum = 0
    # For each dimension, we use a pair of neurons
    for j in range(N):
        # First neuron in the pair has weight +alpha on the j-th input
        neuron_plus_activation = gelu(alpha * x[j])

        # Second neuron has weight -alpha on the j-th input
        neuron_minus_activation = gelu(-alpha * x[j])

        # Add the contribution of this pair to the sum
        hidden_layer_sum += neuron_plus_activation + neuron_minus_activation
    
    # The final output is the weighted sum from the hidden layer
    network_output = output_weight * hidden_layer_sum
    
    return network_output

# --- Demonstration ---
# 1. Define the dimension of the input vector
N = 5
print(f"Input vector dimension (N): {N}")
print(f"Minimum hidden layer width required (H = 2N): {2 * N}")
print("-" * 30)

# 2. Create a random input vector
# Using a small vector to stay in the region where the Taylor approximation is good.
np.random.seed(42)
x_input = np.random.randn(N) * 0.5 
print(f"Input vector x: {x_input}")
print("-" * 30)

# 3. Calculate the true squared norm
y_true = np.linalg.norm(x_input)**2
print(f"True squared norm ||x||^2: {y_true:.6f}")

# 4. Use the network to approximate the squared norm
# A smaller alpha gives a better approximation
alpha = 0.01
y_pred = squared_norm_network(x_input, alpha=alpha)

# Final Equation: network_output = C * sum(g(alpha*x_j) + g(-alpha*x_j))
# C = sqrt(2*pi) / (2 * alpha^2)
C = np.sqrt(2 * np.pi) / (2 * alpha**2)
sum_activations = np.sum([gelu(alpha * x_j) + gelu(-alpha * x_j) for x_j in x_input])

print(f"Network prediction with alpha={alpha}: {y_pred:.6f}")
print("\nCalculation details:")
print(f"Each number in the final equation: y_hat = {C:.4f} * {sum_activations:.10f}")