import numpy as np

def run_variable_step_gd():
    """
    Demonstrates that gradient descent with a variable, summable step size
    can converge to a non-stationary point. This corresponds to a special
    case of algorithm (2).
    """
    # Define a simple linear function f(x) = 2x.
    # The gradient is f'(x) = 2 for all x.
    # A stationary point would require f'(x) = 0, which is impossible.
    f = lambda x: 2 * x
    grad_f = lambda x: 2

    # Initial point
    x_k = 0.0

    # Parameters
    max_iterations = 50
    # A summable step-size sequence
    gamma_k = lambda k: 0.5**(k + 1)

    print(f"Running algorithm (2) on f(x) = 2x with x_0 = 0 and gamma_k = 0.5^(k+1)")
    print("-" * 20)

    total_step = 0
    for k in range(max_iterations):
        step_size = gamma_k(k)
        gradient = grad_f(x_k)
        
        # Update rule for C=R^d: x_{k+1} = x_k - gamma_k * grad_f(x_k)
        x_k_plus_1 = x_k - step_size * gradient
        total_step += step_size * gradient
        
        # print(f"Iteration {k}: x_{k} = {x_k:.6f}, x_{k+1} = {x_k_plus_1:.6f}")
        x_k = x_k_plus_1
        
    # The theoretical sum of the geometric series sum_{k=0 to inf} 0.5^(k+1) is 1.
    # So the theoretical limit point x* is x_0 - grad * sum(gamma_k) = 0 - 2 * 1 = -2.
    limit_point_theory = -2.0
    
    print(f"\nAfter {max_iterations} iterations, the sequence converges to x_k = {x_k:.6f}")
    print(f"The theoretical limit point is x* = {limit_point_theory}")
    
    # Check for stationarity at the limit point
    gradient_at_limit = grad_f(limit_point_theory)
    print(f"The gradient at the limit point x* is f'(x*) = {gradient_at_limit}")

    if gradient_at_limit != 0:
        print("The limit point is not a first-order stationary point.")
    else:
        print("The limit point is a first-order stationary point.")

run_variable_step_gd()