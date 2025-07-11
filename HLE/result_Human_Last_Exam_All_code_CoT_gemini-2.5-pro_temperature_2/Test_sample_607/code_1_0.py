import numpy as np

def f(x):
  """
  A simple smooth function f(x) = x.
  """
  return x

def grad_f(x):
  """
  The gradient of f(x) = x.
  """
  # For a scalar function, gradient is the derivative.
  # For a vector input, this would be a vector of ones.
  return 1.0

def doubly_projected_gradient_descent_demo():
  """
  Demonstrates that Doubly-Projected Gradient Descent can converge
  to a non-stationary point with a summable step-size sequence.

  We use the case where the constraint set C is R, so the projections
  are identity operations. The algorithm becomes x_{k+1} = x_k - gamma_k * grad_f(x_k).
  """
  # Initial point
  x = 0.0
  num_iterations = 50

  print(f"Running Doubly-Projected Gradient Descent for f(x) = x.")
  print(f"Starting at x_0 = {x}")
  print("-" * 30)

  # A summable step-size sequence, gamma_k = (1/2)^(k+1)
  for k in range(num_iterations):
    gamma_k = 0.5**(k + 1)
    gradient = grad_f(x)
    x = x - gamma_k * gradient
    # print(f"k={k}, gamma_k={gamma_k:.4f}, x_{k+1}={x:.6f}") # uncomment for detailed steps

  # The theoretical limit is -1
  x_star = x
  grad_at_star = grad_f(x_star)

  # For an unconstrained problem (C=R^d), the stationarity condition is ||grad f(x*)|| = 0.
  # Proj_{T_{x*}C}(-grad f(x*)) simplifies to -grad f(x*).
  # The norm is |-grad f(x*)|.
  stationarity_value = np.linalg.norm(-grad_at_star)

  print(f"The algorithm converges to the point x* = {x_star:.6f}")
  print(f"The gradient at this point is grad_f(x*) = {grad_at_star}")
  print("\nThe first-order stationarity condition is ||Proj_{T_{x*}C}(-grad f(x*))|| = 0.")
  print("For our unconstrained case, this simplifies to ||-grad f(x*)|| = 0.")
  print("\nFinal equation check:")
  print(f"||-grad f(x*)|| = {stationarity_value}")
  print(f"Is {stationarity_value} equal to 0? {'Yes' if np.isclose(stationarity_value, 0) else 'No'}.")
  print("\nThe algorithm has converged to a non-stationary point.")


if __name__ == '__main__':
  doubly_projected_gradient_descent_demo()
