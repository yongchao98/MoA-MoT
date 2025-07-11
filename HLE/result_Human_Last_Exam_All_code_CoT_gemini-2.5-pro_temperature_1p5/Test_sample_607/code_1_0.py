import numpy as np

def f(x):
    """
    A simple quadratic function f(x) = 0.5 * x^2.
    The gradient is nabla_f(x) = x.
    The only stationary point is x = 0.
    """
    return 0.5 * np.dot(x, x)

def nabla_f(x):
    """Gradient of f(x)."""
    return x

def heavy_ball_method(x_minus_1, x_0, gamma, beta, num_iterations):
    """
    Performs the heavy-ball method iteration.
    x_k+1 = x_k + beta * (x_k - x_k-1) - gamma * nabla_f(x_k)
    """
    print("Heavy-ball method starting from x_0 =", x_0)
    xk = x_0
    xk_minus_1 = x_minus_1
    
    for k in range(num_iterations):
        grad = nabla_f(xk)
        # The core update rule
        xk_plus_1 = xk + beta * (xk - xk_minus_1) - gamma * grad
        
        # Update variables for next iteration
        xk_minus_1 = xk
        xk = xk_plus_1
        
        if k % 10 == 0:
            print(f"Iteration {k}: x = {xk}, f(x) = {f(xk)}")
            
    print("\nFinal point after", num_iterations, "iterations:", xk)
    grad_final = nabla_f(xk)
    print("Gradient at final point:", grad_final)
    if np.linalg.norm(grad_final) < 1e-4:
        print("The algorithm converged to a stationary point.")
    else:
        # Note: For this simple convex function, it will always converge to the stationary point.
        # Demonstrating convergence to a non-stationary point requires a specific non-convex function and parameters.
        print("The algorithm did not converge to a stationary point.")

# Parameters for the simulation
# For the simple quadratic function, these parameters lead to convergence.
initial_point_0 = np.array([10.0, -10.0])
initial_point_minus_1 = np.array([10.1, -10.1])
learning_rate_gamma = 0.1
momentum_beta = 0.9
iterations = 100

# Run the simulation
heavy_ball_method(initial_point_minus_1, initial_point_0, learning_rate_gamma, momentum_beta, iterations)

# (1) x_k+1 = x_k - gamma * nabla_f(x_k)
# (3) x_k+1 = x_k + beta * (x_k - x_k-1) - gamma * nabla_f(x_k)
# Let's show the final update equation for iteration k=1 (the second step)
k = 1
# assume x_0 = [10, -10] and x_-1 = [10.1, -10.1]
x_k_minus_1 = np.array([10.1, -10.1])
x_k = np.array([10.0, -10.0])
grad_k = nabla_f(x_k) # grad_k = [10.0, -10.0]
x_k_plus_1 = x_k + momentum_beta * (x_k - x_k_minus_1) - learning_rate_gamma * grad_k
# x_k+1 = [10.0, -10.0] + 0.9 * ([10.0, -10.0] - [10.1, -10.1]) - 0.1 * [10.0, -10.0]
# x_k+1 = [10.0, -10.0] + 0.9 * [-0.1, 0.1] - [1.0, -1.0]
# x_k+1 = [10.0, -10.0] + [-0.09, 0.09] - [1.0, -1.0]
# x_k+1 = [8.91, -8.91]

print("\nExample of one update step (k=1) for algorithm (3):")
print(f"x_{k+1} = {list(x_k)} + {momentum_beta} * ({list(x_k)} - {list(x_k_minus_1)}) - {learning_rate_gamma} * {list(grad_k)}")
print(f"x_{k+1} = {list(x_k_plus_1)}")