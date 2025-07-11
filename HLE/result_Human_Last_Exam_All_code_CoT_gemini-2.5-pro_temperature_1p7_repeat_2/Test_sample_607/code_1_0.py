import numpy as np

# Define a simple smooth function and its gradient
def f(x):
  """Function f(x) = (x - 10)^2 / 2"""
  return 0.5 * (x - 10)**2

def grad_f(x):
  """Gradient of f(x) which is x - 10"""
  return x - 10

def solve():
    """
    Demonstrates that Gradient Descent converges to a stationary point.
    """
    # Algorithm parameters
    x0 = 0.0      # Initial point
    gamma = 0.1   # Step size
    num_iterations = 100 # Number of iterations

    # Initialize iterate
    x_k = x0

    print(f"Running Gradient Descent for f(x) = (x - 10)^2 / 2")
    print("The stationary point is x = 10, where grad_f(x) = 0.")
    print("-" * 50)
    print(f"Initial point x_0 = {x0:.4f}")
    
    # Run Gradient Descent
    for k in range(num_iterations):
        gradient = grad_f(x_k)
        x_k_plus_1 = x_k - gamma * gradient

        if k == num_iterations - 2:
             # Store values before the last iteration to print the final equation
             final_x_k = x_k
             final_grad = gradient
             final_x_k_plus_1 = x_k_plus_1

        x_k = x_k_plus_1
    
    print(f"Converged point after {num_iterations} iterations: x_k = {x_k:.4f}")
    print(f"Gradient at converged point: grad_f(x_k) = {grad_f(x_k):.4f}")
    print("-" * 50)
    print("The update equation is x_{k+1} = x_k - gamma * grad_f(x_k)")
    print("At the final step, the values are:")
    
    # We output each number in the final equation as requested.
    print(f"{final_x_k_plus_1:.4f} = {final_x_k:.4f} - {gamma} * {final_grad:.4f}")
    print("\nAs the algorithm converges, the gradient approaches 0, and the update step becomes negligible.")

solve()