import math

# Step 1: Define the given constants lambda_1 and lambda_2.
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# Step 2: As per the plan, we assume x_2(0) = 0.
# The expression to be calculated is E = (2/3*lambda1*exp(lambda1/2) - 1/3*exp(lambda1/2))*x_2(0) 
#                                     - 2/3*lambda2*exp(lambda2/2) - 10/3*exp(lambda1/2).
# With x_2(0) = 0, the first term vanishes.
x2_0 = 0

# Step 3: Calculate the terms of the simplified expression.
term1 = (2/3 * lambda1 * math.exp(lambda1/2) - 1/3 * math.exp(lambda1/2)) * x2_0
term2 = -2/3 * lambda2 * math.exp(lambda2/2)
term3 = -10/3 * math.exp(lambda1/2)

# Step 4: Calculate the final value of the expression.
result = term1 + term2 + term3

# Step 5: Print the breakdown of the calculation as requested.
print(f"The expression to calculate is (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2(0) - 2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)")
print(f"Assuming x2(0) = {x2_0}")
print(f"Term 1: (2/3 * {lambda1:.4f} * exp({lambda1/2:.4f}) - 1/3 * exp({lambda1/2:.4f})) * {x2_0:.4f} = {term1:.4f}")
print(f"Term 2: -2/3 * {lambda2:.4f} * exp({lambda2/2:.4f}) = {term2:.4f}")
print(f"Term 3: -10/3 * exp({lambda1/2:.4f}) = {term3:.4f}")
print(f"Final Value = {term1:.4f} + ({term2:.4f}) + ({term3:.4f}) = {result:.4f}")
