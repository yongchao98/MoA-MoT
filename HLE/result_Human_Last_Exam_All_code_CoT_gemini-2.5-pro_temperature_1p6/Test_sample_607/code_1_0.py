import numpy as np

def f(x):
    """The function to minimize: f(x) = (x-10)^2"""
    return (x - 10)**2

def grad_f(x):
    """The gradient of the function: f'(x) = 2(x-10)"""
    return 2 * (x - 10)

def heavy_ball_simulation():
    """
    Simulates the Heavy-ball method on a simple quadratic function.
    """
    # Parameters
    beta = 0.9  # Momentum parameter
    gamma = 0.05  # Step size (learning rate)
    x_minus_1 = 0.0  # Initial point x_{-1}
    x_k = 1.0     # Initial point x_0
    
    # Stationary point
    stationary_point = 10.0

    print("Running Heavy-ball method simulation...")
    print(f"Objective: Minimize f(x) = (x - 10)^2")
    print(f"Stationary point is at x = {stationary_point}")
    print(f"Parameters: beta = {beta}, gamma = {gamma}")
    print(f"Initial points: x_[-1]={x_minus_1}, x_[0]={x_k}\n")
    
    # Store iterates for analysis
    iterates = [x_minus_1, x_k]
    
    for k in range(100):
        # Heavy-ball update rule: x_{k+1} = x_k + beta(x_k - x_{k-1}) - gamma * grad_f(x_k)
        x_k_plus_1 = x_k + beta * (x_k - x_minus_1) - gamma * grad_f(x_k)
        
        # Update points for next iteration
        x_minus_1 = x_k
        x_k = x_k_plus_1
        
        iterates.append(x_k)
        
        # Stop if converged
        if np.abs(x_k - x_minus_1) < 1e-6:
            break
            
    final_x = iterates[-1]
    final_grad = grad_f(final_x)
    
    # For a quadratic, the update is a linear recurrence. We can find the theoretical limit.
    # x_{k+1} = (1 + beta - 2*gamma) * x_k - beta * x_{k-1} + 20*gamma
    # If x_k -> x*, then x* = (1 + beta - 2*gamma)*x* - beta*x* + 20*gamma
    # x* = (1 - 2*gamma)*x* + 20*gamma  => 2*gamma*x* = 20*gamma => x* = 10
    theoretical_limit = (20 * gamma) / (2 * gamma)

    print("Simulation results:")
    print(f"Number of iterations: {k+1}")
    print(f"Final iterate x_k: {final_x}")
    print(f"Stationary point: {stationary_point}")
    print(f"Difference |x_k - stationary_point|: {np.abs(final_x - stationary_point)}")
    print(f"Gradient at final iterate, grad_f(x_k): {final_grad}")
    
    print("\nAs shown by the recurrence relation for this function:")
    print(f"x_k+1 = (1 + {beta} - 2*{gamma})*x_k - {beta}*x_k-1 + 20*{gamma}")
    # We output the numbers of the equation here
    # 1.0 is the coefficient of x* on the left, let's call it c1
    # The coefficient of x* on the right is (1 - 2*gamma)
    c1 = 1.0
    c2 = 1 - 2*gamma
    c3 = 20 * gamma
    print(f"The limit point x* must satisfy the equation: {c1:.4f} * x* = {c2:.4f} * x* + {c3:.4f}")
    print(f"Solving for x* gives the theoretical limit: {theoretical_limit}")
    
heavy_ball_simulation()
<<<H>>>