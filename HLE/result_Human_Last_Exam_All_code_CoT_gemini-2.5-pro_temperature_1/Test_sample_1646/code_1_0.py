import cmath
from scipy.special import gamma
import numpy as np

def f(z):
  """
  Calculates the value of the derived function f(z) = sqrt(π) / Γ(z).
  """
  return cmath.sqrt(cmath.pi) / gamma(z)

def verify_functional_equation(z):
  """
  Verifies that f(z) satisfies the functional equation for a given z.
  f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)
  """
  # Calculate the left-hand side (LHS)
  lhs = f(z)

  # Calculate the components of the right-hand side (RHS)
  term1 = 2**(1 - z)
  term2 = f(z / 2)
  term3 = f((z + 1) / 2)

  # Calculate the right-hand side (RHS)
  rhs = term1 * term2 * term3

  print(f"Verifying the functional equation for z = {z}\n")
  print("The equation is: f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)\n")
  
  # Print each number in the final equation
  print(f"f(z)            = {lhs}")
  print(f"2^(1-z)         = {term1}")
  print(f"f(z/2)          = {term2}")
  print(f"f((z+1)/2)      = {term3}")
  print("-" * 40)
  print(f"LHS: f(z)                          = {lhs}")
  print(f"RHS: 2^(1-z) * f(z/2) * f((z+1)/2) = {rhs}\n")

  # Check if LHS and RHS are close
  if np.isclose(lhs, rhs):
    print("Verification successful: LHS is approximately equal to RHS.")
  else:
    print("Verification failed: LHS is not equal to RHS.")

def verify_initial_condition():
    """
    Verifies that f(1) = sqrt(π).
    """
    val_f1 = f(1)
    val_sqrt_pi = cmath.sqrt(cmath.pi)
    print("\n" + "="*40)
    print("Verifying the initial condition f(1) = sqrt(π)\n")
    print(f"Calculated f(1) = {val_f1}")
    print(f"Expected sqrt(π) = {val_sqrt_pi}\n")
    if np.isclose(val_f1, val_sqrt_pi):
        print("Initial condition satisfied.")
    else:
        print("Initial condition not satisfied.")
    print("="*40)

if __name__ == "__main__":
  # Choose a sample complex number for verification
  z_test = 2.5 + 1.5j
  
  # Verify the functional equation
  verify_functional_equation(z_test)
  
  # Verify the initial condition
  verify_initial_condition()