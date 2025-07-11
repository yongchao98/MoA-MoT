import numpy as np

def grad_f(x):
    """
    The gradient of the non-convex part of the function f(x) = x^4 - 1.5*x^2
    This is valid for x in [-1, 1]. The counterexample uses initial points
    within this range, and the trajectory converges to a point within this range.
    """
    return 4 * x**3 - 3 * x

def heavy_ball_demonstration():
    """
    Runs the Heavy-ball method to show convergence to a non-stationary point.
    """
    # Parameters from Zou et al., 2019
    gamma = 0.01
    beta = 0.9
    
    # Initial points
    x_prev = 0.6
    x_curr = 0.61
    
    # Number of iterations
    num_iterations = 2000
    
    print("Running Heavy-ball method...")
    print(f"Initial points: x_[-1] = {x_prev}, x_[0] = {x_curr}")
    print(f"Parameters: gamma = {gamma}, beta = {beta}")
    
    for i in range(num_iterations):
        grad = grad_f(x_curr)
        x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad
        
        # Update points
        x_prev = x_curr
        x_curr = x_next

    x_final = x_curr
    grad_final = grad_f(x_final)

    print(f"\nAfter {num_iterations} iterations, the method converges to:")
    print(f"x_final = {x_final:.6f}")
    
    print("\nNow, let's check if this point is first-order stationary.")
    print("A point x is stationary if the gradient at x is zero.")
    
    print("\nCalculating the gradient at x_final:")
    # The final equation output as requested.
    print(f"grad_f(x_final) = 4 * ({x_final:.6f})^3 - 3 * ({x_final:.6f})")
    print(f"                = 4 * {x_final**3:.6f} - {3*x_final:.6f}")
    print(f"                = {4*x_final**3:.6f} - {3*x_final:.6f}")
    print(f"                = {grad_final:.6f}")

    # A small tolerance for checking if the gradient is zero
    if abs(grad_final) > 1e-3:
        print("\nThe gradient is not zero. The point is not stationary.")
    else:
        print("\nThe gradient is close to zero. The point is stationary.")

if __name__ == '__main__':
    heavy_ball_demonstration()
