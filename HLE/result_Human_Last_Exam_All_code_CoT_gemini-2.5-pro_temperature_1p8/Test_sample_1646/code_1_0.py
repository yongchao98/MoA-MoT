import numpy as np
from scipy.special import gamma

def f(z):
  """
  The derived function f(z) = sqrt(pi) / Gamma(z).
  """
  return np.sqrt(np.pi) / gamma(z)

def verify_solution():
  """
  Verifies the solution by testing the functional equation and the initial condition.
  """
  print("--- Verification of the functional equation f(z) = 2^(1-z) * f(z/2) * f((z+1)/2) ---")
  
  # Use a sample complex number for z
  z = 2.5 + 1.5j
  
  # Calculate the Left-Hand Side (LHS)
  lhs = f(z)
  
  # Calculate the Right-Hand Side (RHS)
  term_f_z2 = f(z / 2)
  term_f_z12 = f((z + 1) / 2)
  term_pow = 2**(1 - z)
  rhs = term_pow * term_f_z2 * term_f_z12
  
  print(f"For z = {z}:")
  print(f"LHS = f({z}) = {lhs}")
  print(f"RHS = 2^(1-({z})) * f({z/2}) * f({(z+1)/2})")
  print(f"RHS = ({term_pow}) * ({term_f_z2}) * ({term_f_z12})")
  print(f"RHS = {rhs}")
  
  # Check if LHS and RHS are close (accounting for floating-point inaccuracies)
  is_close = np.allclose(lhs, rhs)
  print(f"\nAre LHS and RHS approximately equal? {is_close}")
  
  print("\n" + "-"*60 + "\n")
  
  print("--- Verification of the initial condition f(1) = sqrt(pi) ---")
  
  # Calculate f(1)
  val_f1 = f(1)
  val_sqrt_pi = np.sqrt(np.pi)
  
  print(f"Calculated value: f(1) = {val_f1}")
  print(f"Expected value: sqrt(pi) = {val_sqrt_pi}")
  
  is_close_cond = np.allclose(val_f1, val_sqrt_pi)
  print(f"\nAre the calculated and expected values approximately equal? {is_close_cond}")

  print("\n" + "-"*60 + "\n")
  print("The explicit form of the function is: f(z) = sqrt(pi) / Gamma(z)")

if __name__ == '__main__':
  verify_solution()
