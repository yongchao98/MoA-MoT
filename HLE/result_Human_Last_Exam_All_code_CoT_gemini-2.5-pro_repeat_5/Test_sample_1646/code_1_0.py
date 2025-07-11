import cmath
# We use scipy's gamma function as it is a standard, reliable implementation.
from scipy.special import gamma as gamma_func

# Based on the derivation, the explicit form of the function is f(z) = sqrt(pi) / Gamma(z).
# This script will numerically verify this solution.

def f(z):
  """
  The derived function f(z) = sqrt(pi) / Gamma(z).
  """
  # scipy.special.gamma can handle complex numbers.
  return cmath.sqrt(cmath.pi) / gamma_func(z)

# --- Numerical Verification ---
print("--- Numerical Verification ---")

# 1. Verify the functional equation: f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)
z = 1.5 + 2j # A sample complex number for verification
lhs = f(z)
rhs = (2**(1 - z)) * f(z / 2) * f((z + 1) / 2)

print(f"Verifying the functional equation for z = {z}:")
print(f"LHS = f(z) = {lhs}")
print(f"RHS = 2^(1-z) * f(z/2) * f((z+1)/2) = {rhs}")
print(f"Are LHS and RHS approximately equal? {cmath.isclose(lhs, rhs)}\n")

# 2. Verify the initial condition: f(1) = sqrt(pi)
val_at_1 = f(1)
expected_val_at_1 = cmath.sqrt(cmath.pi)

print("Verifying the initial condition f(1) = sqrt(pi):")
print(f"Calculated f(1) = {val_at_1}")
print(f"Expected sqrt(pi) = {expected_val_at_1}")
print(f"Are they approximately equal? {cmath.isclose(val_at_1, expected_val_at_1)}\n")

# --- Final Answer ---
print("--- Final Answer ---")
print("The explicit form of the function is: f(z) = sqrt(pi) / Gamma(z)")
print("The numerical constant in this final equation is pi.")
# As requested, outputting the number(s) in the final equation.
# The primary constant here is pi.
print("Value of the number pi:")
print(cmath.pi)