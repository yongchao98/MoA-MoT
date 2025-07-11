import numpy as np

def f(x):
    """
    A simple non-convex function: f(x) = (x^2 - 4)^2 = x^4 - 8x^2 + 16
    Stationary points are at x=0, x=2, x=-2.
    """
    return (x**2 - 4)**2

def grad_f(x):
    """
    Gradient of f(x)
    """
    return 4 * x * (x**2 - 4)

def heavy_ball_method(start_x, gamma, beta, max_iter=1000, tol=1e-6):
    """
    Implements the Heavy-ball method.
    """
    # Initialization
    x_prev = float(start_x)
    x_curr = float(start_x) # Start with zero momentum
    
    print(f"Starting Heavy-ball method at x = {start_x}")
    print(f"Parameters: gamma = {gamma}, beta = {beta}")
    
    for i in range(max_iter):
        x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad_f(x_curr)
        
        # Check for convergence
        if abs(x_next - x_curr) < tol:
            print(f"\nConverged after {i+1} iterations.")
            x_curr = x_next
            break
            
        x_prev = x_curr
        x_curr = x_next
    else:
        print("\nDid not converge within max iterations.")

    # Final point and gradient
    final_x = x_curr
    final_grad = grad_f(final_x)
    
    print(f"Final point x^* = {final_x}")
    print(f"Gradient at final point grad_f(x^*) = {final_grad}")
    
    # Check if the limit point is stationary
    is_stationary = np.isclose(final_grad, 0)
    print(f"Is the limit point stationary? {is_stationary}")
    if is_stationary:
        print("As demonstrated, the algorithm converged to a stationary point.")
    else:
        print("The algorithm converged to a NON-STATIONARY point.")


# Run the demonstration
# Parameters are chosen for stable convergence
heavy_ball_method(start_x=1.0, gamma=0.01, beta=0.9)
