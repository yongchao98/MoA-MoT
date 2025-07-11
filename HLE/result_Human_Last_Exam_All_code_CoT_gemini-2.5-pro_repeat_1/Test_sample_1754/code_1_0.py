import math

# Define the constants lambda1 and lambda2
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# The expression to calculate is:
# (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2(0) 
# - 2/3 * lambda2 * exp(lambda2/2) 
# - 10/3 * exp(lambda1/2)

# We assume x2(0) = 0, which is standard for reachability problems.
x2_0 = 0

# Calculate each term of the expression.
# The first term becomes 0 because it is multiplied by x2(0).
term1 = (2/3 * lambda1 * math.exp(lambda1 / 2) - 1/3 * math.exp(lambda1 / 2)) * x2_0
term2 = - (2/3) * lambda2 * math.exp(lambda2 / 2)
term3 = - (10/3) * math.exp(lambda1 / 2)

# The final value is the sum of these terms.
final_value = term1 + term2 + term3

# As requested, output each number in the final equation.
# The equation is: final_value = term1 + term2 + term3
print("Assuming x2(0) = 0, the final equation is the sum of the following terms:")
print(f"Term 1: {term1}")
print(f"Term 2: {term2}")
print(f"Term 3: {term3}")
print("\nFinal Equation:")
print(f"{term1} + ({term2}) + ({term3}) = {final_value}")

# Print the final numerical answer.
print(f"\nThe final calculated value is: {final_value}")