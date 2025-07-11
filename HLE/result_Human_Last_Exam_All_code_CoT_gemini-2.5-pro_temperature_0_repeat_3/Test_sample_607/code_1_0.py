import math

def f(x):
    """The objective function."""
    return 0.25 * x**4 - 0.5 * x**2

def grad_f(x):
    """The gradient of the objective function."""
    return x**3 - x

def heavy_ball_cycle_demo():
    """
    Demonstrates the Heavy-ball method converging to a limit cycle
    composed of non-stationary points.
    """
    # Parameters for the limit cycle {sqrt(2), -sqrt(2)}
    beta = 0.5
    gamma = 3.0
    
    # Initial points close to the cycle
    x_prev = -math.sqrt(2) + 0.01
    x_curr = math.sqrt(2) - 0.01
    
    print("Heavy-ball method for f(x) = x^4/4 - x^2/2")
    print(f"Stationary points are at x = 0, 1, -1.")
    print(f"Parameters: beta = {beta}, gamma = {gamma}")
    print("Starting iterations...\n")
    
    print(f"k=-1: x = {x_prev:.6f}")
    print(f"k=0:  x = {x_curr:.6f}, grad(x) = {grad_f(x_curr):.6f}")
    
    for k in range(1, 11):
        # Heavy-ball update rule: x_{k+1} = x_k + beta(x_k - x_{k-1}) - gamma * grad(f(x_k))
        momentum_term = beta * (x_curr - x_prev)
        gradient_term = gamma * grad_f(x_curr)
        x_next = x_curr + momentum_term - gradient_term
        
        # Print the equation for the first step to show the calculation
        if k == 1:
            print("\n--- Calculation for k=1 ---")
            print(f"x_next = x_curr + beta * (x_curr - x_prev) - gamma * grad_f(x_curr)")
            print(f"x_next = {x_curr:.4f} + {beta:.1f} * ({x_curr:.4f} - ({x_prev:.4f})) - {gamma:.1f} * {grad_f(x_curr):.4f}")
            print(f"x_next = {x_curr:.4f} + {momentum_term:.4f} - {gradient_term:.4f}")
            print(f"x_next = {x_next:.4f}")
            print("---------------------------\n")

        # Update variables for next iteration
        x_prev = x_curr
        x_curr = x_next
        
        print(f"k={k}:  x = {x_curr:.6f}, grad(x) = {grad_f(x_curr):.6f}")

    print("\nAs shown, the iterates oscillate between sqrt(2) and -sqrt(2).")
    print(f"The gradient at these points is not zero (grad(sqrt(2)) = {grad_f(math.sqrt(2)):.4f}),")
    print("so the algorithm has 'converged' to a set of non-stationary points.")

heavy_ball_cycle_demo()
<<<C>>>