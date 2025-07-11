import math

# Step 1: Define the given constants lambda_1 and lambda_2.
sqrt_17 = math.sqrt(17)
lambda_1 = (1 + sqrt_17) / 2
lambda_2 = (1 - sqrt_17) / 2

print(f"Given lambda_1 = (1 + sqrt(17))/2, its value is: {lambda_1}")
print(f"Given lambda_2 = (1 - sqrt(17))/2, its value is: {lambda_2}")
print("-" * 20)

# The expression to evaluate is:
# (2/3 * lambda_1 * exp(lambda_1/2) - 1/3 * exp(lambda_1/2)) * x_2(0) - 2/3 * lambda_2 * exp(lambda_2/2) - 10/3 * exp(lambda_1/2)

# Step 2: Apply the assumption that the initial state is zero.
# In controllability problems, we often analyze reachability from the origin, meaning x(0) = 0.
# Thus, x_2(0) = 0.
x_2_0 = 0
print(f"Assuming the initial state x(0) is zero, we have x_2(0) = {x_2_0}")
print("The expression simplifies because the first term becomes zero.")
print("-" * 20)

# Step 3: Calculate the simplified expression.
# The simplified expression is: - (2/3) * lambda_2 * exp(lambda_2/2) - (10/3) * exp(lambda_1/2)

# Calculate each part of the equation
exp_lambda_1_half = math.exp(lambda_1 / 2)
exp_lambda_2_half = math.exp(lambda_2 / 2)

term1 = -(2/3) * lambda_2 * exp_lambda_2_half
term2 = -(10/3) * exp_lambda_1_half

final_value = term1 + term2

# Step 4: Print the components of the final equation and the result.
print("The final equation is: result = (-2/3 * lambda_2 * exp(lambda_2/2)) + (-10/3 * exp(lambda_1/2))")
print(f"Value of lambda_2: {lambda_2}")
print(f"Value of exp(lambda_2 / 2): {exp_lambda_2_half}")
print(f"Value of the first term (-2/3 * lambda_2 * exp(lambda_2/2)): {term1}")
print("-" * 20)
print(f"Value of lambda_1: {lambda_1}")
print(f"Value of exp(lambda_1 / 2): {exp_lambda_1_half}")
print(f"Value of the second term (-10/3 * exp(lambda_1/2)): {term2}")
print("-" * 20)
print(f"The final calculated value is the sum of these two terms.")
print(f"Final Value = {term1} + {term2} = {final_value}")