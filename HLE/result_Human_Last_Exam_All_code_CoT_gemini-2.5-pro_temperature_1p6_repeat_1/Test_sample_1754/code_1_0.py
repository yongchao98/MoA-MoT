import math

# Step 1: Define the constants lambda_1 and lambda_2
sqrt_17 = math.sqrt(17)
lambda1 = (1 + sqrt_17) / 2
lambda2 = (1 - sqrt_17) / 2

# Step 2: As stated in the plan, we assume x_2(0) = 0.
# The expression simplifies to - (2/3)*lambda2*e^(lambda2/2) - (10/3)*e^(lambda1/2).
x2_0 = 0

# Step 3: Calculate each term of the full expression for the output string.
term1_coeff_factor1 = (2/3) * lambda1 * math.exp(lambda1 / 2)
term1_coeff_factor2 = -(1/3) * math.exp(lambda1 / 2)
term1_coeff = term1_coeff_factor1 + term1_coeff_factor2
term1 = term1_coeff * x2_0

term2 = -(2/3) * lambda2 * math.exp(lambda2 / 2)
term3 = -(10/3) * math.exp(lambda1 / 2)

# Step 4: Calculate the final value
final_value = term1 + term2 + term3

# Step 5: Print the equation with the numbers and the final result.
# Note: The problem asks to output each number in the final equation.
# Full expression: ( (2/3)*lambda1*e^(lambda1/2) - (1/3)*e^(lambda1/2) ) * x2_0 - (2/3)*lambda2*e^(lambda2/2) - (10/3)*e^(lambda1/2)
print("Assuming the system starts from the origin, we have x_2(0) = 0.")
print("The expression to calculate is:")
print(f"( (2/3) * {lambda1:.4f} * e^({lambda1:.4f}/2) - (1/3) * e^({lambda1:.4f}/2) ) * {x2_0} - (2/3) * ({lambda2:.4f}) * e^({lambda2:.4f}/2) - (10/3) * e^({lambda1:.4f}/2)")
print("\nPlugging in the values for the exponential terms:")
print(f"( {term1_coeff:.4f} ) * {x2_0} + ( {term2:.4f} ) + ( {term3:.4f} )")
print(f"\nSimplified calculation: {term1:.4f} + {term2:.4f} + {term3:.4f}")
print(f"\nFinal Value = {final_value:.4f}")