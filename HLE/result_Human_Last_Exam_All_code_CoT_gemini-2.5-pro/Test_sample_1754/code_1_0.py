import math

# Step 1: Define the constants lambda_1 and lambda_2
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# The expression to calculate is:
# (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2(0) - (2/3) * lambda2 * exp(lambda2/2) - (10/3) * exp(lambda1/2)
# We assume this is a reachability problem from the origin, so x(0) = 0, which implies x2(0) = 0.
# The expression simplifies to: - (2/3) * lambda2 * exp(lambda2/2) - (10/3) * exp(lambda1/2)

# Step 2: Calculate each term of the simplified expression
term1 = -(2/3) * lambda2 * math.exp(lambda2 / 2)
term2 = -(10/3) * math.exp(lambda1 / 2)

# Step 3: Output each number (term) in the final equation
print(f"Assuming x2(0) = 0, the expression simplifies to the sum of two terms.")
print(f"First term: -(2/3) * lambda2 * exp(lambda2/2) = {term1}")
print(f"Second term: -(10/3) * exp(lambda1/2) = {term2}")

# Step 4: Calculate the final value
final_value = term1 + term2

# Print the final result
print("\nThe final value is:")
print(final_value)