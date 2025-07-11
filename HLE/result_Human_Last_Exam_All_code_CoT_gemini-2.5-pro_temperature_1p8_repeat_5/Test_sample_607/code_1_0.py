import numpy as np

# This script demonstrates that algorithm (2), Doubly-projected gradient descent,
# can converge to a point that is not first-order stationary.

# Let's consider a simple 1D case.
# Function: f(x) = x. This is a smooth function.
# Gradient: nabla_f(x) = 1. This is never zero, so there are no stationary points.
def grad_f(x):
  """Computes the gradient of f(x) = x."""
  return 1.0

# Constraint Set: C = R (the entire real line).
# In this case, Proj_C is the identity function.
def proj_C(y):
  """Projection onto the set C=R."""
  return y

# Tangent Cone: T_x C = R for any x in C.
# So, Proj_{T_x C} is also the identity function.
def proj_tangent_cone(v, x):
  """Projection onto the tangent cone T_x C for C=R."""
  return v

# Algorithm (2) update rule:
# x_{k+1} = Proj_C(x_k + gamma_k * Proj_{T_{x_k}C}(-nabla f(x_k)))

# With our choices, this simplifies to:
# x_{k+1} = x_k + gamma_k * (-1) = x_k - gamma_k

# We choose a summable step-size sequence, gamma_k = 0.5^(k+1), so the
# algorithm converges. The infinite sum of these step sizes is 1.0.
x_k = 0.0  # Initial point x_0
num_iterations = 50

print("Demonstrating Algorithm (2) for f(x) = x and C = R.")
print(f"The gradient is always 1, so no stationary points exist.")
print(f"Starting at x_0 = {x_k}")

# Run the iteration
for k in range(num_iterations):
  gamma_k = 0.5**(k + 1)
  neg_grad = -grad_f(x_k)
  proj_neg_grad = proj_tangent_cone(neg_grad, x_k)
  x_k = proj_C(x_k + gamma_k * proj_neg_grad)

# The theoretical limit is x_0 - sum(gamma_k) = 0 - 1 = -1.
limit_point = x_k
print(f"\nAfter {num_iterations} iterations, the algorithm converges to x = {limit_point}")

# Check the stationarity condition at the limit point.
# A point x* is stationary if ||Proj_{T_{x*}C}(-nabla f(x*))|| = 0.
# For our case, this is ||-nabla f(x*)|| = 0, which means nabla f(x*) = 0.
grad_at_limit = grad_f(limit_point)
print(f"The gradient at this limit point is: {grad_at_limit}")
stationarity_norm = np.linalg.norm(proj_tangent_cone(-grad_at_limit, limit_point))
print(f"The norm of the projected gradient is: {stationarity_norm}")
print("\nConclusion: The algorithm has converged to a point, but this point is not stationary.")