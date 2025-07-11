import math

# Step 1: Define the given constants lambda_1 and lambda_2.
lambda_1 = (1 + math.sqrt(17)) / 2
lambda_2 = (1 - math.sqrt(17)) / 2

# Step 2: Assume x_2(0) = 0 as the initial state is not specified.
# This simplifies the expression to be calculated.
x2_0 = 0

# Step 3: Define the components of the expression.
term1_coeff_val = (2/3 * lambda_1 * math.exp(lambda_1 / 2) - 1/3 * math.exp(lambda_1 / 2))
term1 = term1_coeff_val * x2_0
term2 = (2/3) * lambda_2 * math.exp(lambda_2 / 2)
term3 = (10/3) * math.exp(lambda_1 / 2)

# Step 4: Calculate the final value.
# The expression is term1 - term2 - term3
final_value = term1 - term2 - term3

# Step 5: Print the equation with the calculated values.
print(f"The expression to calculate is:")
print(f"({term1_coeff_val}) * {x2_0} - ({term2}) - ({term3})")
print(f"This simplifies to:")
print(f"{term1} - {term2} - {term3} = {final_value}")
