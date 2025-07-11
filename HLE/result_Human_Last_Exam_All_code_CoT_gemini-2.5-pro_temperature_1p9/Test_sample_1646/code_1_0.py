import cmath
import math
from scipy.special import gamma

def f(z):
  """
  Calculates the value of the derived function f(z) = sqrt(pi) / Gamma(z).
  """
  return cmath.sqrt(math.pi) / gamma(z)

print("The determined explicit form of the function is f(z) = sqrt(π) / Γ(z).")
print("We will now verify that this function satisfies the given conditions.")
print("-" * 60)

# 1. Verify the initial condition f(1) = sqrt(π)
z1 = 1
f_at_1 = f(z1)
sqrt_pi_val = cmath.sqrt(math.pi)

print(f"Step 1: Verifying the initial condition f(1) = sqrt(π):")
print(f"  The number in the condition is 1.")
print(f"  Calculated f({z1}) = {f_at_1}")
print(f"  Expected sqrt(π) = {sqrt_pi_val}")
print(f"  Is the condition met? {cmath.isclose(f_at_1, sqrt_pi_val)}")
print("-" * 60)

# 2. Verify the functional equation for a sample complex value, z = 3.5 + 2j
z_val = 3.5 + 2j

# Calculate the Left Hand Side (LHS): f(z)
lhs = f(z_val)

# Calculate the Right Hand Side (RHS): 2^(1 - z) * f(z/2) * f((z+1)/2)
# The numbers in the equation are 2, 1
base = 2
exponent = 1 - z_val
coeff = base ** exponent
term_f1 = f(z_val / base)
term_f2 = f((z_val + 1) / base)
rhs = coeff * term_f1 * term_f2

print(f"Step 2: Verifying the functional equation for z = {z_val}:")
print(f"  The functional equation is f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)")
print(f"  The numbers in the equation are 2 and 1.")

print("\n  Calculating Left Hand Side (LHS):")
print(f"    LHS = f({z_val})")
print(f"    LHS = {lhs}")

print("\n  Calculating Right Hand Side (RHS):")
print(f"    RHS = ({base}^({exponent})) * f({z_val / base}) * f({(z_val + 1) / base})")
print(f"    RHS = {rhs}")

print("\n  Comparison:")
print(f"    LHS and RHS are close: {cmath.isclose(lhs, rhs)}")