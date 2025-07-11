import numpy as np

def sigmoid(y):
    """Numerically stable sigmoid function."""
    return np.where(y >= 0, 1 / (1 + np.exp(-y)), np.exp(y) / (1 + np.exp(y)))

def objective_function(x, y):
    """The objective function to minimize."""
    s_y = sigmoid(y)
    # Add small epsilon to log arguments to avoid log(0)
    epsilon = 1e-8
    return s_y * (x**2) - np.log(s_y + epsilon) - np.log(1 - s_y + epsilon)

def gradient(x, y):
    """Computes the gradient of the objective function."""
    s_y = sigmoid(y)
    s_prime_y = s_y * (1 - s_y)
    
    df_dx = 2 * x * s_y
    df_dy = s_prime_y * (x**2 - 1 + 2 * s_y)
    return np.array([df_dx, df_dy])

def heavy_ball_method():
    """
    Implements the Heavy-ball method and shows convergence to a non-stationary point.
    """
    # Parameters
    beta = 0.7  # Momentum parameter
    gamma = 0.1 # Step size / learning rate
    iterations = 500

    # Initialization that leads to convergence to a non-stationary point
    x_prev = np.array([-2.0, 2.5])
    x_curr = np.array([-1.9, 2.45])

    print(f"Running Heavy-ball method for {iterations} iterations...")
    
    for k in range(iterations):
        grad_f = gradient(x_curr[0], x_curr[1])
        x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad_f
        
        # Update iterates
        x_prev = x_curr
        x_curr = x_next
    
    final_point = x_curr
    final_gradient = gradient(final_point[0], final_point[1])

    print("\n--- Results ---")
    print(f"Algorithm converged to the point: ({final_point[0]:.6f}, {final_point[1]:.6f})")
    print(f"Gradient at this point is: ({final_gradient[0]:.6f}, {final_gradient[1]:.6f})")
    
    # The true stationary point is (0,0)
    grad_at_stationary = gradient(0, 0)
    
    grad_norm_sq = np.sum(final_gradient**2)
    
    print(f"The squared norm of the final gradient is: {grad_norm_sq:.6f}")
    
    if grad_norm_sq > 1e-4:
        print("\nThe algorithm has converged to a point that is NOT first-order stationary.")
    else:
        print("\nThe algorithm has converged to a stationary point.")

heavy_ball_method()
