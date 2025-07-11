import numpy as np
from scipy.special import gamma

def f(z):
  """
  Calculates the value of the function f(z) = sqrt(pi) / Gamma(z).
  """
  return np.sqrt(np.pi) / gamma(z)

def verify_functional_equation(z):
  """
  Verifies the functional equation f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)
  for a given complex number z.
  """
  print(f"--- Verifying for z = {z} ---")

  # Calculate the Left-Hand Side (LHS)
  lhs = f(z)

  # Calculate each component of the Right-Hand Side (RHS)
  term1 = 2**(1 - z)
  term2 = f(z / 2)
  term3 = f((z + 1) / 2)
  
  # Calculate the full RHS
  rhs = term1 * term2 * term3
  
  print("The functional equation is: f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)\n")

  print("--- Numerical values of each term ---")
  print(f"f(z)           = {lhs}")
  print(f"2^(1-z)        = {term1}")
  print(f"f(z/2)         = {term2}")
  print(f"f((z+1)/2)     = {term3}")
  print("-" * 20)
  print(f"LHS f(z)             = {lhs}")
  print(f"RHS (product of terms) = {rhs}")
  
  # Check if LHS and RHS are close (allowing for floating point inaccuracies)
  if np.isclose(lhs, rhs):
    print("\nConclusion: The functional equation holds true.")
  else:
    print("\nConclusion: The functional equation does not hold.")

def verify_initial_condition():
    """
    Verifies the initial condition f(1) = sqrt(pi).
    """
    print("\n--- Verifying initial condition f(1) = sqrt(pi) ---")
    val_f1 = f(1)
    val_sqrt_pi = np.sqrt(np.pi)
    print(f"Calculated f(1) = {val_f1}")
    print(f"Expected sqrt(pi) = {val_sqrt_pi}")
    if np.isclose(val_f1, val_sqrt_pi):
        print("\nConclusion: The initial condition holds true.")
    else:
        print("\nConclusion: The initial condition does not hold.")

if __name__ == "__main__":
  # Choose a sample complex number for z
  z_test = 2.5 + 1.5j
  
  # Perform the verifications
  verify_functional_equation(z_test)
  verify_initial_condition()

  print("\n\nThe explicit form of the function is f(z) = sqrt(pi) / Gamma(z)")
