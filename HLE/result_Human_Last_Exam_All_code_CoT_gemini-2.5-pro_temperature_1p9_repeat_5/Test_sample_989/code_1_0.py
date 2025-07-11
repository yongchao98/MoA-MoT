import numpy as np
from scipy.optimize import minimize
import warnings

# Suppress overflow warnings for the demonstration
warnings.filterwarnings("ignore", message="overflow encountered in exp")
warnings.filterwarnings("ignore", message="Values in x were outside bounds during a minimize step")

def f(x):
  """
  A strictly convex function f(x) = e^x.
  Its second derivative is e^x, which is always > 0, so it's strictly convex.
  """
  return np.exp(x[0])

def demonstrate_no_minimum():
  """
  Demonstrates that a strictly convex function may not have a global minimum.
  We attempt to use a numerical optimizer to find the minimum of f(x) = e^x.
  The optimizer will fail because the function's value continues to decrease
  as x approaches negative infinity.
  """
  # Starting point for the optimization
  x0 = [0]
  
  # Attempt to minimize the function
  # The method 'BFGS' is a standard gradient-based optimization algorithm.
  # For a convex function, if it converges, it finds the global minimum.
  result = minimize(f, x0, method='BFGS')
  
  print("--- Optimization Attempt for f(x) = e^x ---")
  print(f"Statement to test: 'Any strictly convex function has a unique global minimizer.'")
  print(f"Counterexample function: f(x) = e^x, which is strictly convex.")
  print("\nResult from scipy.optimize.minimize:")
  print(f"Success: {result.success}")
  print(f"Status Message: {result.message}")
  print(f"Final function value (fun): {result.fun}")
  print(f"Final x value (x): {result.x[0]}")
  
  print("\n--- Conclusion ---")
  if not result.success:
    print("The optimization failed to converge to a minimum.")
    print("This demonstrates that although f(x) = e^x is strictly convex, it does not have a global minimum.")
    print("Therefore, the statement is FALSE.")
  else:
    # This case is unlikely but included for completeness.
    print("The optimizer claimed success, but the result should be inspected.")
    print("For f(x)=e^x, any result is non-sensical as the true infimum is 0 at x = -infinity.")

# Run the demonstration
demonstrate_no_minimum()
