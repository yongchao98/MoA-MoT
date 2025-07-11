import numpy as np
from scipy.special import gamma

# The derived explicit form of the function is f(z) = sqrt(pi) / Γ(z).
# This script will verify this solution.

# In the final equation, f(z) = C / Γ(z), the constant C is sqrt(pi).
# We will print the values of these constants.
pi_val = np.pi
C = np.sqrt(pi_val)

print("--- The Final Equation's Components ---")
print(f"The explicit form is f(z) = C / Γ(z), where Γ is the Gamma function.")
print(f"The value of the number pi is: {pi_val}")
print(f"The value of the constant C = sqrt(pi) is: {C}")
print("-" * 40)

def f(z):
    """
    The derived function f(z) = sqrt(pi) / Gamma(z).
    """
    return np.sqrt(np.pi) / gamma(z)

# --- Verification ---
print("\n--- Verifying the Solution ---\n")

# 1. Verify the functional equation: f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)
# We choose a random complex number to test this.
z_test = 3.5 + 2.5j

# Calculate the left-hand side (LHS) of the equation
lhs = f(z_test)

# Calculate the right-hand side (RHS) of the equation
rhs = (2**(1 - z_test)) * f(z_test / 2) * f((z_test + 1) / 2)

print(f"1. Testing the functional equation for z = {z_test}:")
print(f"   LHS = f(z)                      = {lhs}")
print(f"   RHS = 2^(1-z)f(z/2)f((z+1)/2)   = {rhs}")
print(f"   Are the two sides equal? {np.isclose(lhs, rhs)}\n")


# 2. Verify the initial condition: f(1) = sqrt(pi)
f_at_1 = f(1)
sqrt_pi = np.sqrt(np.pi)

print("2. Testing the condition at z = 1:")
print(f"   Calculated f(1) = {f_at_1}")
print(f"   Expected sqrt(pi) = {sqrt_pi}")
print(f"   Are the two values equal? {np.isclose(f_at_1, sqrt_pi)}\n")

print("-" * 40)
print("\nThe verifications confirm that the derived function is correct.")
