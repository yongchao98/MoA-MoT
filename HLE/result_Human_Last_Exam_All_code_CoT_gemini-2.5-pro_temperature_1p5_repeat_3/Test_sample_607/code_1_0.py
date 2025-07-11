import numpy as np

def f(x):
  """The objective function."""
  return 0.5 * x**2

def grad_f(x):
  """The gradient of the objective function."""
  return x

def run_counterexample():
    """
    This function provides a counterexample for Algorithm (2), showing it can
    converge to a non-stationary point.

    We use f(x) = 0.5*x^2 and C = R. The only stationary point is x=0.
    The algorithm simplifies to x_{k+1} = x_k - gamma_k * x_k.
    We choose a summable step-size sequence gamma_k = 1/2^(k+2).
    """
    x0 = 10.0  # Initial point
    num_iterations = 50  # Sufficient for convergence

    print("--- Counterexample for Algorithm (2) ---")
    print(f"Objective function: f(x) = 0.5 * x^2")
    print(f"Gradient: grad_f(x) = x")
    print("The only stationary point is x = 0.\n")
    
    x = x0
    
    # Demonstrate the first step of the iteration
    k=0
    gamma_k = 1.0 / (2**(k + 2))
    grad_k = grad_f(x)
    x_next = x - gamma_k * grad_k
    
    print("Equation: x_{k+1} = x_k - gamma_k * grad_f(x_k)")
    print("First iteration (k=0):")
    print(f"x_0 = {x}")
    print(f"gamma_0 = 1/2^(0+2) = {gamma_k}")
    print(f"grad_f(x_0) = {grad_k}")
    print(f"x_1 = {x:.4f} - {gamma_k:.4f} * {grad_k:.4f} = {x_next:.4f}\n")
    
    # Run the algorithm for many iterations
    for i in range(num_iterations):
      gamma_i = 1.0 / (2**(i + 2))
      x = x * (1 - gamma_i)

    x_star = x
    grad_f_x_star = grad_f(x_star)

    print(f"Starting at x_0 = {x0:.4f}")
    print(f"After {num_iterations} iterations, the sequence converges to x_star = {x_star:.4f}")
    print(f"The gradient at this limit point is grad_f(x_star) = {grad_f_x_star:.4f}")

    # Check if the limit point is stationary
    if not np.isclose(grad_f_x_star, 0.0):
      print("\nThe limit point x_star is NOT a stationary point because its gradient is not zero.")
    else:
      print("\nThe limit point is a stationary point.")

run_counterexample()