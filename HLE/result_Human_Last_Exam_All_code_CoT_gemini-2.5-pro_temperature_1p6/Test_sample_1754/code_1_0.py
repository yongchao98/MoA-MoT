import math

# Define the constants lambda_1 and lambda_2
sqrt_17 = math.sqrt(17)
lambda_1 = (1 + sqrt_17) / 2
lambda_2 = (1 - sqrt_17) / 2

# The expression to calculate with x2(0)=0
x2_0 = 0
term1_coeff_val = (2/3) * lambda_1 * math.exp(lambda_1 / 2) - (1/3) * math.exp(lambda_1 / 2)
term2_val = -(2/3) * lambda_2 * math.exp(lambda_2 / 2)
term3_val = -(10/3) * math.exp(lambda_1 / 2)

# Calculate the final value
final_value = term1_coeff_val * x2_0 + term2_val + term3_val

# Print the final equation with all numbers substituted
print(f"({term1_coeff_val:.4f}) * {x2_0} + ({term2_val:.4f}) + ({term3_val:.4f}) = {final_value:.4f}")
print("The final equation is:")
print(f"({(2/3)*lambda_1:.4f}*exp({lambda_1/2:.4f}) - {(1/3)*math.exp(lambda_1/2):.4f})*{x2_0} - {(2/3)*lambda_2:.4f}*exp({lambda_2/2:.4f}) - {(10/3)*math.exp(lambda_1/2):.4f}")
print("Which simplifies to:")
print(f"0 - {(2/3)*lambda_2:.4f}*exp({lambda_2/2:.4f}) - {(10/3)*math.exp(lambda_1/2):.4f} = {final_value:.4f}")

# Just print the final value as requested for the final answer block
# Do not use this for the explanation part
final_answer = final_value