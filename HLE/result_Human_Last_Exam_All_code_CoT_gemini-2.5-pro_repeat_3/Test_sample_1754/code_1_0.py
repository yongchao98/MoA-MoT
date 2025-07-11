import math

# Step 1: Define the given eigenvalues lambda_1 and lambda_2.
lambda_1 = (1 + math.sqrt(17)) / 2
lambda_2 = (1 - math.sqrt(17)) / 2

# Step 2: State the assumption for x_2(0).
# The problem is underdetermined as x_2(0) is not specified.
# A common approach in controllability is to assume the system starts from the origin.
# We will assume x_1(0) = 0 and x_2(0) = 0.
x_2_0 = 0

# Step 3: Calculate each term in the expression based on the assumption.
# The expression is:
# (2/3 * lambda_1 * exp(lambda_1/2) - 1/3 * exp(lambda_1/2)) * x_2(0)
# - (2/3) * lambda_2 * exp(lambda_2/2)
# - (10/3) * exp(lambda_1/2)

# With x_2(0) = 0, the first term is 0.
term1 = (2/3 * lambda_1 * math.exp(lambda_1 / 2) - 1/3 * math.exp(lambda_1 / 2)) * x_2_0
term2 = - (2/3) * lambda_2 * math.exp(lambda_2 / 2)
term3 = - (10/3) * math.exp(lambda_1 / 2)

# Step 4: Calculate the final result.
result = term1 + term2 + term3

# Step 5: Print the final equation with the calculated values of each term.
print("Based on the assumption that the system starts from the origin, x_2(0) = 0.")
print("The expression to be calculated is:")
print(f"({(2/3 * lambda_1 * math.exp(lambda_1 / 2) - 1/3 * math.exp(lambda_1 / 2)):.4f}) * {x_2_0} - (2/3 * {lambda_2:.4f} * exp({lambda_2 / 2:.4f})) - (10/3 * exp({lambda_1 / 2:.4f}))")
print("The final equation with numerical values is:")
print(f"{term1:.4f} + ({term2:.4f}) + ({term3:.4f}) = {result:.4f}")
print("\nFinal Result:")
print(result)
