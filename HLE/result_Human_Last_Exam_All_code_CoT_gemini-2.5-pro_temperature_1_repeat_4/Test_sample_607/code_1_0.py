import numpy as np

def f(x, y):
    """The objective function."""
    return 0.5 * (x - 1)**2 - 2 * (y - x**2)**2

def grad_f(x, y):
    """The gradient of the objective function."""
    df_dx = (x - 1) - 8 * x * (y - x**2)
    df_dy = -4 * (y - x**2)
    return np.array([df_dx, df_dy])

def heavy_ball_example():
    """
    Demonstrates the Heavy-ball method potentially converging to a non-stationary point.
    This is for illustrative purposes; convergence can be sensitive to parameters and initialization.
    """
    # Parameters from literature suggesting convergence to non-stationary point (0,0)
    # The actual convergence is a delicate matter.
    beta = 1.0
    gamma = 0.045 
    
    # Initialization near the non-stationary point (0,0)
    # The point (0,0) lies on the parabola y=x^2
    x_prev = np.array([0.1, 0.1**2])
    x_curr = np.array([0.05, 0.05**2])
    
    print("Heavy-ball method example for f(x,y) = 1/2*(x-1)^2 - 2*(y-x^2)^2")
    print(f"Targeting non-stationary point (0,0). Gradient at (0,0) is {grad_f(0,0)}.")
    print(f"Parameters: beta = {beta}, gamma = {gamma}")
    print(f"Initial x_k-1: {x_prev}")
    print(f"Initial x_k:   {x_curr}")
    print("-" * 20)

    for k in range(20):
        grad = grad_f(x_curr[0], x_curr[1])
        x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad
        
        # Update iterates
        x_prev = x_curr
        x_curr = x_next
        
        if (k + 1) % 5 == 0:
            print(f"Iteration {k+1}: x = {x_curr}, grad_norm = {np.linalg.norm(grad):.4f}")

    print("-" * 20)
    print(f"Final position after 20 iterations: {x_curr}")
    print("For specific non-convex functions and parameters, the Heavy-ball method (3) has been shown to be able to converge to a non-stationary point.")
    print("Under standard assumptions (e.g., constant step size, or diminishing step sizes with an infinite sum), methods (1) and (2) are guaranteed to have limit points that are stationary.")
    print("Therefore, only method (3) stands out as having this property under standard conditions.")

heavy_ball_example()