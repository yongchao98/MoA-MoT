import math

# Define the eigenvalues based on the provided formulas
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# Based on the steady-state analysis, x2(0) = 0.
x2_0 = 0

# The expression to evaluate is:
# (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2_0 - 2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)
# With x2_0 = 0, the expression simplifies to:
# -2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)

# Define the coefficients and terms of the final equation
c1 = -2/3
c2 = -10/3
t_val = 1/2

# Calculate the terms
term1 = c1 * lambda2 * math.exp(lambda2 * t_val)
term2 = c2 * math.exp(lambda1 * t_val)

print("Based on a steady-state analysis, we found x2(0) = 0.")
print("The expression to evaluate simplifies to: c1 * lambda2 * exp(lambda2 * t) + c2 * exp(lambda1 * t)")
print("\nIndividual numbers in the final equation:")
print(f"c1 = {c1}")
print(f"lambda2 = {lambda2}")
print(f"t = {t_val}")
print(f"c2 = {c2}")
print(f"lambda1 = {lambda1}")

print(f"\nCalculating the two terms:")
print(f"First term (c1 * lambda2 * exp(lambda2 * t)): {term1}")
print(f"Second term (c2 * exp(lambda1 * t)): {term2}")

# Calculate the final result
result = term1 + term2

print(f"\nThe final result is the sum of these two terms:")
print(f"Result = {result}")
