import math

# Define the constants lambda1 and lambda2
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# Assumption for the initial state x2(0)
x2_0 = 0

# Calculate each part of the expression
# First part: coefficient of x2(0)
term1_coeff = (2/3 * lambda1 * math.exp(lambda1 / 2)) - (1/3 * math.exp(lambda1 / 2))

# Second part: term with lambda2
term2 = (2/3 * lambda2 * math.exp(lambda2 / 2))

# Third part: term with lambda1
term3 = (10/3 * math.exp(lambda1 / 2))

# Calculate the final result
result = term1_coeff * x2_0 - term2 - term3

# Print the numbers in the final equation
print("The equation to solve is:")
print("(a) * x2(0) - (b) - (c)")
print(f"where a = {term1_coeff}")
print(f"      x2(0) = {x2_0} (by assumption)")
print(f"      b = {term2}")
print(f"      c = {term3}")
print("\nFinal Equation:")
print(f"({term1_coeff}) * ({x2_0}) - ({term2}) - ({term3})")
print(f"\nResult: {result}")