import numpy as np
from scipy.special import gamma

# The constant in our derived equation is C = sqrt(pi)
C = np.sqrt(np.pi)

def f(z):
  """
  This is the derived explicit form for f(z), which is
  f(z) = sqrt(pi) / Gamma(z).
  """
  return C / gamma(z)

# --- Verification Step 1: Initial Condition ---
print("Verifying the initial condition: f(1) = sqrt(pi)")
f_of_1 = f(1)
sqrt_pi = np.sqrt(np.pi)
print(f"Calculated f(1) = {f_of_1}")
print(f"Expected sqrt(pi) = {sqrt_pi}")
print(f"Is f(1) equal to sqrt(pi)? {np.isclose(f_of_1, sqrt_pi)}\n")

# --- Verification Step 2: Functional Equation ---
# We use a sample complex number to test the equation
z = 1.5 + 2j
print(f"Verifying the functional equation for z = {z}")

# Calculate the Left-Hand Side (LHS)
lhs = f(z)
# Calculate the Right-Hand Side (RHS)
rhs = (2**(1 - z)) * f(z / 2) * f((z + 1) / 2)

print(f"LHS of the equation, f(z) = {lhs}")
print(f"RHS of the equation, 2^(1-z) * f(z/2) * f((z+1)/2) = {rhs}")
print(f"Are LHS and RHS approximately equal? {np.isclose(lhs, rhs)}\n")

# --- Final Answer: The Explicit Form ---
print("The explicit form of the function is determined to be:")
# The final equation is f(z) = C / Gamma(z)
# Printing the numbers/constants in the final equation.
print(f"f(z) = C / Gamma(z), where C is the square root of pi.")
print(f"C = {C}")
# Printing the full formula again
print(f"\nFinal Equation: f(z) = sqrt({np.pi}) / Gamma(z)")