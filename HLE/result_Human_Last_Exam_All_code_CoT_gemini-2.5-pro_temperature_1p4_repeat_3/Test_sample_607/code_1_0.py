import numpy as np

def f(x):
    """The function to minimize: f(x) = 0.5*x^2 + x."""
    return 0.5 * x**2 + x

def grad_f(x):
    """The gradient of the function: f'(x) = x + 1."""
    return x + 1

def demonstrate_convergence_to_non_stationary_point():
    """
    Demonstrates that a gradient-based method with a specific choice of
    diminishing step sizes can converge to a non-stationary point.
    """
    # The algorithm is a special case of (2) with C=R^d:
    # x_{k+1} = x_k - gamma_k * grad_f(x_k)
    
    # Target point: x_star = 0.
    # Note that grad_f(0) = 1, so x_star = 0 is NOT a stationary point.
    
    # We choose a sequence of iterates x_k = 1/k and derive the
    # required step sizes gamma_k to produce this sequence.
    # From x_{k+1} = x_k - gamma_k * (x_k + 1), we can solve for gamma_k:
    # gamma_k = (x_k - x_{k+1}) / (x_k + 1)
    # Substituting x_k = 1/k, we get gamma_k = 1/((k+1)*(k+2)).
    
    # Number of iterations
    num_iterations = 100
    
    # Initial point, chosen as x_1 = 1/1 = 1
    x = 1.0
    
    print(f"Targeting convergence to non-stationary point x=0.")
    print(f"Stationary point is at x=-1, where grad_f(x)=0.")
    print(f"Gradient at target point x=0 is grad_f(0) = {grad_f(0)}")
    print("-" * 50)
    print(f"Iter 1: x = {x:.6f}, grad_f(x) = {grad_f(x):.6f}")

    for k in range(1, num_iterations + 1):
        # Calculate the special step size for this iteration
        # gamma_k = 1 / ((k+1)*(k+2))
        gamma = 1.0 / ((k + 1.0) * (k + 2.0))
        
        # Calculate the gradient at the current point
        gradient = grad_f(x)
        
        # Apply the update rule
        x = x - gamma * gradient
        
        # Print progress at certain iterations
        if k+1 in [2, 3, 4, 5, 10, 50, 100]:
             print(f"Iter {k+1}: x = {x:.6f}, grad_f(x) = {grad_f(x):.6f}")
             
    print("-" * 50)
    print(f"Final state after {num_iterations+1} iterations:")
    print(f"x = {x:.6f}")
    print(f"Analytical limit for x_k is 0.")
    print(f"Gradient at final x: grad_f(x) = {grad_f(x):.6f}")
    print(f"The gradient converges to 1, not 0.")

demonstrate_convergence_to_non_stationary_point()