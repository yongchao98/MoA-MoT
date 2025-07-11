import math

# Define the numbers in the final expression
n1 = 10000
n2 = 2
n3 = 2
n4 = 7

# Calculate the dominant root lambda
lambda_val = n3 + math.sqrt(n4)

# Calculate the limit L
limit_L = n2 * math.log(lambda_val)

# Calculate the final value
final_value = math.floor(n1 * limit_L)

# The problem asks to output the numbers in the final equation.
# The final equation is floor(10000 * 2 * ln(2 + sqrt(7))).
print(f"The expression to compute is floor({n1} * {n2} * log({n3} + sqrt({n4}))).")
print(f"The final integer result is: {final_value}")