import numpy as np

def heavy_ball_optimizer():
    """
    This function demonstrates the Heavy-ball optimization algorithm.
    It minimizes the function f(x) = (x - 10)^2.
    The stationary point is at x=10, where the gradient is 0.
    """
    
    # Parameters for the Heavy-ball method
    beta = 0.9  # Momentum parameter
    gamma = 0.05 # Step size (learning rate)
    
    # Initial points
    x_prev = 0.0
    x_curr = 0.0
    
    # Objective function and its gradient
    f = lambda x: (x - 10)**2
    grad_f = lambda x: 2 * (x - 10)
    
    print("Running Heavy-ball method to minimize f(x) = (x - 10)^2")
    print(f"Initial points: x_[-1] = {x_prev}, x_[0] = {x_curr}\n")
    
    # Iteration loop
    num_iterations = 50
    for k in range(num_iterations):
        # Calculate the gradient at the current point
        gradient = grad_f(x_curr)
        
        # Calculate the momentum term
        momentum_term = beta * (x_curr - x_prev)
        
        # Apply the Heavy-ball update rule
        x_next = x_curr + momentum_term - gamma * gradient
        
        # Update points for the next iteration
        x_prev = x_curr
        x_curr = x_next
        
    print(f"After {num_iterations} iterations:")
    print(f"The algorithm converged to x = {x_curr:.6f}")
    
    final_gradient = grad_f(x_curr)
    print(f"The gradient at this point is nabla f(x) = {final_gradient:.6f}")
    
    print("\nIn this typical case for a convex function, the algorithm converges to the stationary point x=10 where the gradient is 0.")
    print("However, for some non-convex functions, it is claimed that this method can converge to a non-stationary point.")

    # Retrieve the values from the last iteration to display the final equation
    # We need to re-calculate the final momentum term and gradient used
    final_momentum = beta * (x_curr - x_prev)
    final_gradient_val = grad_f(x_prev) # grad was calculated at x_curr (which is now x_prev)
    
    print("\nThe final update step equation was:")
    print(f"x_{num_iterations} = x_{num_iterations-1} + beta * (x_{num_iterations-1} - x_{num_iterations-2}) - gamma * grad_f(x_{num_iterations-1})")
    print(f"{x_curr:.4f} = {x_prev:.4f} + {beta:.1f} * ({x_prev:.4f} - x_prev_before_last(x_prev, final_momentum, beta):.4f}) - {gamma:.2f} * {final_gradient_val:.4f}")

def x_prev_before_last(x_prev, momentum, beta):
    # Helper to reconstruct x_{k-2} from x_{k-1} and momentum_k
    if beta == 0: return "N/A"
    return x_prev - momentum / beta

heavy_ball_optimizer()