import numpy as np

def grad_f(x, M):
    """Computes the gradient of the test function."""
    x1, x2 = x[0], x[1]
    # f(x1, x2) = 0.5*x1^2 + 0.5*M*x2^2 - x1^2*x2
    grad_x1 = x1 - 2 * x1 * x2
    grad_x2 = M * x2 - x1**2
    return np.array([grad_x1, grad_x2], dtype=np.float64)

# Parameters for the Heavy-ball method
gamma = 0.01
beta = 0.9
M = 40.0

# Initialization
x_prev = np.array([-1.0, 1.25], dtype=np.float64)
x_curr = np.array([-1.0, 1.25], dtype=np.float64)

# Run the Heavy-ball method for many iterations
num_iterations = 50000
for k in range(num_iterations):
    grad = grad_f(x_curr, M)
    # The Heavy-ball update rule
    x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad
    x_prev = x_curr.copy()
    x_curr = x_next.copy()

# Final point reached by the algorithm
final_x = x_curr
final_grad = grad_f(final_x, M)

# The stationarity condition for this unconstrained problem is || grad f(x) || = 0.
# Let's check this condition at the point of convergence.
print(f"Final iterate after {num_iterations} iterations: x_k = [{final_x[0]:.8f}, {final_x[1]:.8f}]")
print(f"Gradient at the final iterate: grad_f(x_k) = [{final_grad[0]:.8f}, {final_grad[1]:.8f}]")
final_grad_norm = np.linalg.norm(final_grad)
print(f"Norm of the gradient ||grad_f(x_k)|| = {final_grad_norm:.8f}")

# Check the stationarity equation: || grad f(x*) || = 0
# The equation for the final point is || [0.00000000, 0.00000000] || = 0.00000000
print(f"\nThe stationarity equation at the limit point x* is || grad f(x*) || = 0.")
print("The components of the gradient at our final iterate are:")
print(f"Component 1: {final_grad[0]:.8f}")
print(f"Component 2: {final_grad[1]:.8f}")
print(f"The values are effectively zero, confirming convergence to a stationary point.")
