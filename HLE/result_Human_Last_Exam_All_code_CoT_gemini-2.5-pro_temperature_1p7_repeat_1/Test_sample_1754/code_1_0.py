import math

# Step 1: Define the constants lambda_1 and lambda_2
sqrt_17 = math.sqrt(17)
lambda_1 = (1 + sqrt_17) / 2
lambda_2 = (1 - sqrt_17) / 2
T = 0.5 # Corresponds to the term 1/2 in the exponent

# Step 2: The expression to evaluate is:
# ( (2/3)*lambda_1*exp(lambda_1/2) - (1/3)*exp(lambda_1/2) )*x_2(0) 
# - (2/3)*lambda_2*exp(lambda_2/2) - (10/3)*exp(lambda_1/2)
# Based on the problem's formulation as a controllability problem, we assume x(0) = 0, so x_2(0) = 0.
# The expression simplifies to: - (2/3)*lambda_2*exp(lambda_2/2) - (10/3)*exp(lambda_1/2)

x2_0 = 0

# Step 3: Calculate each numerical part of the original equation
term1_factor1_coeff1 = 2/3
term1_factor1_lambda1 = lambda_1
term1_factor1_exp = math.exp(lambda_1 * T)
term1_factor2_coeff1 = -1/3
term1_factor2_exp = math.exp(lambda_1 * T)
term1_factor1 = term1_factor1_coeff1 * term1_factor1_lambda1 * term1_factor1_exp
term1_factor2 = term1_factor2_coeff1 * term1_factor2_exp
term1 = term1_factor1 + term1_factor2

term2_coeff = -2/3
term2_lambda2 = lambda_2
term2_exp = math.exp(lambda_2 * T)
term2 = term2_coeff * term2_lambda2 * term2_exp

term3_coeff = -10/3
term3_exp = math.exp(lambda_1 * T)
term3 = term3_coeff * term3_exp

# Step 4: Calculate the final result
result = term1 * x2_0 + term2 + term3

# Step 5: Print the final equation with each number substituted
print("Assuming x_2(0) = 0, the equation becomes:")
print(f"( ({term1_factor1_coeff1:.4f} * {term1_factor1_lambda1:.4f} * {term1_factor1_exp:.4f}) + ({term1_factor2_coeff1:.4f} * {term1_factor2_exp:.4f}) ) * {x2_0} + ({term2_coeff:.4f} * {term2_lambda2:.4f} * {term2_exp:.4f}) + ({term3_coeff:.4f} * {term3_exp:.4f})")
print("\nWhich simplifies to:")
print(f"0 + ({term2:.4f}) + ({term3:.4f})")
print(f"\nFinal Result: {result:.4f}")
