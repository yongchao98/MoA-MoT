import numpy as np

def run_heavy_ball_counterexample():
    """
    This function implements the Heavy-ball method for the counter-example 
    proposed by Zavriev and Kostyuk (1993).

    The literature claims that for this function and specific parameters, the 
    algorithm can converge to a non-stationary point.
    """

    # 1. Define the function and its gradient
    # The function is f(x) = (8/3)x^3 - 4x^2 + 2x
    # The gradient is f'(x) = 8x^2 - 8x + 2 = 2 * (2x - 1)^2
    def grad_f(x):
        return 2 * (2 * x - 1)**2

    # The stationary point is where the gradient is zero:
    # 2 * (2x - 1)^2 = 0  => x = 0.5
    stationary_point = 0.5
    
    # The claim is that the method can converge to x = 1, which is not stationary.
    non_stationary_target = 1.0
    grad_at_target = grad_f(non_stationary_target)
    
    print(f"Stationary point is at x = {stationary_point}")
    print(f"The non-stationary point of interest is x = {non_stationary_target}")
    print(f"The gradient at x = {non_stationary_target} is {grad_at_target:.2f}, which is not zero.")
    print("-" * 30)

    # 2. Set Heavy-ball parameters and initial points
    # These are the parameters and initial points cited in relation to the counter-example.
    beta = 3.0 / 4.0
    gamma = 1.0 / 2.0
    
    # Initial points
    x_prev = 0.0  # x_{k-1}
    x_curr = 1.0  # x_k

    print("Running Heavy-ball iterations...")
    print(f"Parameters: beta = {beta}, gamma = {gamma}")
    print(f"Initial points: x_(-1) = {x_prev}, x_0 = {x_curr}\n")
    print(f"k = -1, x = {x_prev:.6f}")
    print(f"k =  0, x = {x_curr:.6f}")

    # 3. Run the iteration
    for k in range(1, 11):
        # Calculate the gradient at the current point
        gradient = grad_f(x_curr)
        
        # Heavy-ball update rule: 
        # x_{k+1} = x_k + beta * (x_k - x_{k-1}) - gamma * grad_f(x_k)
        x_next = x_curr + beta * (x_curr - x_prev) - gamma * gradient

        # Prepare for the next iteration
        x_prev = x_curr
        x_curr = x_next
        
        print(f"k = {k:2d}, x = {x_curr:.6f}")

    print("\nNote: For the given initial conditions, the sequence appears to diverge,")
    print("not converge to 1. However, the theory states that there exist")
    print("initial conditions for which convergence to a non-stationary point occurs.")


run_heavy_ball_counterexample()