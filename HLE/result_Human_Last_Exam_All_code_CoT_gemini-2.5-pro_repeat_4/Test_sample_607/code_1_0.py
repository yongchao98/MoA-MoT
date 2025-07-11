import numpy as np

def f(x):
    """The function to minimize."""
    return 0.5 * (x - 10)**2

def grad_f(x):
    """The gradient of the function."""
    return x - 10

def demonstrate_convergence_to_non_stationary():
    """
    Demonstrates that algorithm (2) can converge to a non-stationary point.
    For this demonstration, we use C = R^d, so the algorithm simplifies to
    gradient descent with a variable step size gamma_k.
    """
    # Initial point
    x_k = 0.0
    
    # Stationary point
    stationary_point = 10.0
    
    print(f"The function is f(x) = 0.5 * (x - 10)^2")
    print(f"The unique stationary point is x = {stationary_point}")
    print(f"Starting at x_0 = {x_k}")
    print("-" * 30)

    # Number of iterations
    num_iterations = 50
    
    for k in range(num_iterations):
        # Summable step size
        gamma_k = 0.5**(k + 1)
        
        # Calculate gradient
        gradient = grad_f(x_k)
        
        # Update rule
        x_k_plus_1 = x_k - gamma_k * gradient
        
        # Uncomment the line below to see each step
        # print(f"k={k}, x_k={x_k:.6f}, gamma_k={gamma_k:.6f}, grad={gradient:.6f}, x_{k+1}={x_k_plus_1:.6f}")

        # Check for convergence
        if np.abs(x_k_plus_1 - x_k) < 1e-9:
            print(f"Converged after {k+1} iterations.")
            x_k = x_k_plus_1
            break
            
        x_k = x_k_plus_1

    print(f"\nThe algorithm converged to the point x* = {x_k:.8f}")
    
    # Check if the limit point is stationary
    final_gradient = grad_f(x_k)
    print(f"The gradient at this point is grad_f(x*) = {final_gradient:.8f}")
    
    if np.abs(final_gradient) > 1e-6:
        print("\nThe limit point is NOT a stationary point.")
    else:
        print("\nThe limit point IS a stationary point.")

demonstrate_convergence_to_non_stationary()