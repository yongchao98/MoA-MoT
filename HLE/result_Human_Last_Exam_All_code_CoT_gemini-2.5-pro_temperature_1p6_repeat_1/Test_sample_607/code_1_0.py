import numpy as np

def heavy_ball_method_example():
    """
    Demonstrates the Heavy-ball method on a simple function.
    """
    # Objective function f(x) = 0.5 * x[0]**2 + 5 * x[1]**2
    # The gradient is grad(f)(x) = [x[0], 10*x[1]]
    # The unique stationary point (and minimizer) is [0, 0].
    
    def grad_f(x):
        return np.array([x[0], 10 * x[1]])

    # Algorithm parameters
    gamma = 0.05  # Learning rate
    beta = 0.9    # Momentum parameter
    
    # Initialization
    x_prev = np.array([10.0, 2.0]) # x_{k-1}
    x_curr = np.array([9.5, 1.5])  # x_k
    
    iterations = 50

    print("Running Heavy-ball method...")
    print(f"Initial x0 = {x_prev}")
    print(f"Initial x1 = {x_curr}")
    
    for k in range(1, iterations + 1):
        # Calculate gradient at the current point
        grad = grad_f(x_curr)
        
        # Calculate momentum term
        momentum_term = beta * (x_curr - x_prev)
        
        # Calculate gradient step
        gradient_step = gamma * grad
        
        # Update the position
        x_next = x_curr + momentum_term - gradient_step
        
        # Update previous and current points for the next iteration
        x_prev = x_curr
        x_curr = x_next
        
        if k % 10 == 0:
            print(f"Iteration {k}: x = {x_curr}, grad = {grad_f(x_curr)}")

    print("\n--- Final Result ---")
    print(f"Point reached after {iterations} iterations: x_{iterations+1} = {x_curr}")
    
    final_grad = grad_f(x_prev)
    final_momentum = beta * (x_prev - x_curr) # this is approx -beta*(x_curr-x_prev)
    
    # The last update was x_curr = x_prev + beta*(x_prev - x_before_prev) - gamma*grad(x_prev)
    # Let's show the final calculated update equation
    # x_{k+1} = x_k + beta * (x_k - x_{k-1}) - gamma * nabla_f(x_k)
    # The final point x_curr was x_{iterations+1}. It was calculated from x_{iterations} and x_{iterations-1}.
    # We don't have x_before_prev saved easily, so let's just print the values for the last computed step.
    
    print("Final update equation terms:")
    print(f"x_{iterations+1} = x_{iterations} + beta * (x_{iterations} - x_{iterations-1}) - gamma * grad(f(x_{iterations}))")
    print(f"{x_curr} = {x_prev} + {beta} * ({x_prev} - [previous x]) - {gamma} * {final_grad}")
    print(f"Note: In this standard case, as the method converges, the point approaches the stationary point [0, 0] where the gradient is [0, 0].")


heavy_ball_method_example()