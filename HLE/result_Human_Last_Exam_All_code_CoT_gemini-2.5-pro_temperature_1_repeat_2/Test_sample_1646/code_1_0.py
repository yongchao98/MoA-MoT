import numpy as np
from scipy.special import gamma

def f(z):
  """
  Calculates the value of the derived function f(z) = sqrt(pi) / Gamma(z).
  """
  return np.sqrt(np.pi) / gamma(z)

# --- Verification ---

# 1. Verify the initial condition f(1) = sqrt(pi)
val_at_1 = f(1)
sqrt_pi_val = np.sqrt(np.pi)

print("--- Verification of the initial condition f(1) = sqrt(pi) ---")
print(f"Calculated f(1): {val_at_1}")
print(f"Value of sqrt(pi): {sqrt_pi_val}")
print(f"Is f(1) == sqrt(pi)? {np.isclose(val_at_1, sqrt_pi_val)}\n")

# 2. Verify the functional equation for a sample complex number.
# The functional equation is: f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2)
# We will test this with z = 2.5 + 1.5j.
z = 2.5 + 1.5j

# Calculate the Left-Hand Side (LHS) of the equation
lhs = f(z)

# Calculate the Right-Hand Side (RHS) of the equation
# The numbers in the equation's structure are 2 and 1.
coeff_base = 2
coeff_exp_const = 1
rhs_factor1 = coeff_base**(coeff_exp_const - z)
rhs_factor2 = f(z / 2)
rhs_factor3 = f((z + 1) / 2)
rhs = rhs_factor1 * rhs_factor2 * rhs_factor3

print(f"--- Verification of the functional equation for z = {z} ---")
print(f"Equation: f(z) = 2^(1 - z) * f(z/2) * f((z+1)/2)\n")

print("Left-Hand Side (LHS):")
print(f"f({z}) = {lhs}\n")

print("Right-Hand Side (RHS):")
print("The components of the RHS are:")
print(f"  The number '2' raised to the power of (1 - z), i.e., 2^(1-{z}), is: {rhs_factor1}")
print(f"  The term f(z/2) for z/2 = {z/2} is: {rhs_factor2}")
print(f"  The term f((z+1)/2) for (z+1)/2 = {(z+1)/2} is: {rhs_factor3}")
print(f"Total RHS value = {rhs}\n")

print("--- Comparison ---")
print(f"LHS = {lhs}")
print(f"RHS = {rhs}")
print(f"Are LHS and RHS approximately equal? {np.isclose(lhs, rhs)}")