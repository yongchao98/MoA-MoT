import numpy as np

def f(x):
    """A simple linear function."""
    return x

def grad_f(x):
    """The gradient of the function f(x)."""
    # The gradient is constant, regardless of x.
    return 1.0

def run_algorithm_2_counterexample(x0, num_iterations):
    """
    Runs a simplified version of Algorithm (2) that demonstrates
    convergence to a non-stationary point.
    
    In R^1 with C=R, the algorithm simplifies to x_{k+1} = x_k - gamma_k * grad_f(x_k).
    """
    x = x0
    print(f"Initial point x_0 = {x}")
    
    for k in range(num_iterations):
        # Choose a step size that decays to zero
        gamma_k = 0.5**(k + 1)
        
        # Calculate the gradient
        gradient = grad_f(x)
        
        # Update step
        x = x - gamma_k * gradient
        
        if k < 10 or (k+1) % 10 == 0:
            print(f"Iteration {k+1:3d}: gamma_k = {gamma_k:.6f}, x_{k+1} = {x:.6f}")
            
    return x

def main():
    """
    Main function to run the demonstration.
    """
    initial_x = 1.0
    iterations = 50
    
    print("Demonstrating that Algorithm (2) can converge to a non-stationary point.")
    print("We use f(x) = x, so grad_f(x) = 1 for all x.")
    print("The algorithm converges to x* = 0, but grad_f(0) = 1, so it's not a stationary point.")
    print("-" * 50)
    
    final_x = run_algorithm_2_counterexample(initial_x, iterations)
    
    print("-" * 50)
    
    # The theoretical limit point
    x_star = 0.0
    
    # The gradient at the limit point
    grad_at_x_star = grad_f(x_star)
    
    # For C=R^d, the stationarity condition is ||grad_f(x*)|| = 0
    stationarity_value = np.linalg.norm(grad_at_x_star)
    
    print(f"The simulation ran for {iterations} iterations.")
    print(f"The final iterate x_k is: {final_x:.8f}")
    print(f"The theoretical limit point x* is: {x_star}")
    print("\nChecking the stationarity condition at the limit point x*:")
    # The "final equation" is the stationarity condition
    print(f"The gradient at the limit point is: grad_f(x*) = {grad_at_x_star}")
    print(f"The stationarity condition is ||grad_f(x*)|| = 0.")
    print(f"Our calculated value is: ||{grad_at_x_star}|| = {stationarity_value}")
    
    if not np.isclose(stationarity_value, 0):
        print("\nConclusion: The limit point is NOT a stationary point.")
    else:
        print("\nConclusion: The limit point is a stationary point.")

if __name__ == "__main__":
    main()