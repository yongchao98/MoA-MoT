import math

# Step 1: Define the constants lambda_1 and lambda_2
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# Step 2: Assume the initial state x2(0) is 0 based on the problem's context
x2_0 = 0

# Step 3: Calculate the numerical value of each part of the expression
# The expression is: (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2(0) - 2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)

# Coefficient of x2(0)
coeff_x2_0 = (2/3) * lambda1 * math.exp(lambda1 / 2) - (1/3) * math.exp(lambda1 / 2)

# The term with x2(0)
term1 = coeff_x2_0 * x2_0

# The second term
term2_val = (2/3) * lambda2 * math.exp(lambda2 / 2)

# The third term
term3_val = (10/3) * math.exp(lambda1 / 2)

# Step 4: Calculate the final result
final_value = term1 - term2_val - term3_val

# Step 5: Print the equation with each number plugged in
print("The equation with the calculated values is:")
print(f"({coeff_x2_0}) * {x2_0} - ({term2_val}) - ({term3_val}) = {final_value}")
print("\nSimplified calculation:")
print(f"0 - {term2_val} - {term3_val} = {final_value}")
