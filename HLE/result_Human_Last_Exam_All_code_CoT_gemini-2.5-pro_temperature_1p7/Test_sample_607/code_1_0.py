import numpy as np

def heavy_ball_demonstration():
    """
    Demonstrates that for the Heavy-ball method to converge to a point,
    the gradient at that point must be zero.
    """
    # Define a simple smooth function f(x) = 0.5 * (x - 5)^2
    # The unique stationary point is x = 5.
    f = lambda x: 0.5 * (x - 5)**2
    grad_f = lambda x: x - 5
    stationary_point = 5

    # Heavy-ball parameters
    beta = 0.9  # Momentum parameter
    gamma = 0.05  # Learning rate
    
    # Initialization
    x_prev = 0.0
    x_curr = 0.1
    
    print("Running Heavy-ball method for f(x) = 0.5 * (x - 5)^2")
    print("Stationary point is at x = 5.")
    print("-" * 60)
    print(f"{'Iteration':<10} {'x_k':<15} {'grad_f(x_k)':<15} {'update_residual':<15}")
    print("-" * 60)

    for k in range(100):
        # Calculate gradient
        gradient = grad_f(x_curr)
        
        # Heavy-ball update
        x_next = x_curr + beta * (x_curr - x_prev) - gamma * gradient
        
        # My analysis showed that gamma * grad_f(x_k) should be equal to
        # x_k + beta*(x_k - x_{k-1}) - x_{k+1}.
        # Let's call the right side the "update residual".
        # Note: here we compute it using x_k, x_{k-1} and x_{k+1}.
        update_residual = x_curr + beta * (x_curr - x_prev) - x_next
        
        if k % 10 == 0:
            print(f"{k:<10} {x_curr:<15.6f} {gradient:<15.6f} {update_residual:<15.6f}")

        # Update iterates for next loop
        x_prev = x_curr
        x_curr = x_next
    
    print("-" * 60)
    print(f"Final iterate x_k: {x_curr:.6f}")
    print(f"Gradient at final iterate grad_f(x_k): {grad_f(x_curr):.6f}")
    
    print("\nAs the iterates x_k converge to the stationary point 5,")
    print("both the gradient and the 'update_residual' converge to 0.")
    print("For a sequence to converge to a limit x*, the update residual")
    print("x_k + beta*(x_k - x_{k-1}) - x_{k+1} must go to 0.")
    print("This forces gamma*grad_f(x_k) to go to 0, which implies that the limit point")
    print("must be a stationary point.")

heavy_ball_demonstration()

# Final Answer based on the step-by-step analysis
# None of the algorithms can converge to a point that is not first-order stationary.
final_answer = "H"
print(f"\nFinal Answer: {final_answer}")