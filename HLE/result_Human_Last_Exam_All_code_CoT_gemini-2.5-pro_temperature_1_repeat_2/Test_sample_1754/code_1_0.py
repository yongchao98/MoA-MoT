import math

# Step 1: Define the constants lambda_1 and lambda_2
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# Step 2: Assume x2(0) = 0, as is common in reachability problems
x2_0 = 0

# Step 3: Calculate the terms of the expression
term1_coeff_part1 = (2/3) * lambda1 * math.exp(lambda1 / 2)
term1_coeff_part2 = (-1/3) * math.exp(lambda1 / 2)
term1_coeff = term1_coeff_part1 + term1_coeff_part2

term1 = term1_coeff * x2_0

term2_coeff = (-2/3) * lambda2
term2 = term2_coeff * math.exp(lambda2 / 2)

term3_coeff = (-10/3)
term3 = term3_coeff * math.exp(lambda1 / 2)

# Step 4: Calculate the final value
final_value = term1 + term2 + term3

# Print the equation with the calculated values
print(f"The expression to calculate is:")
print(f"({term1_coeff_part1:.4f} + {term1_coeff_part2:.4f}) * {x2_0} + ({term2_coeff:.4f} * exp({lambda2/2:.4f})) + ({term3_coeff:.4f} * exp({lambda1/2:.4f}))")
print(f"Which simplifies to:")
print(f"{term1:.4f} + {term2:.4f} + {term3:.4f} = {final_value:.4f}")

# Print the final value in the required format
print("\nFinal Value:")
print(final_value)