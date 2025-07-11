import math

# Step 1: Define the given constants lambda_1 and lambda_2
sqrt_17 = math.sqrt(17)
lambda1 = (1 + sqrt_17) / 2
lambda2 = (1 - sqrt_17) / 2

# Step 2: State the assumption for x2(0)
# Based on the analysis, for the problem to have a unique solution,
# we assume the system starts from rest, which means x2(0) = 0.
x2_0 = 0

# Step 3: Define the expression to calculate E.
# E = ( (2/3)*l1*exp(l1/2) - (1/3)*exp(l1/2) ) * x2(0) 
#     - (2/3)*l2*exp(l2/2) 
#     - (10/3)*exp(l1/2)
# With x2(0) = 0, the expression simplifies.

# Step 4: Calculate each term of the expression as requested.
# The instruction is: "Remember in the final code you still need to output each number in the final equation!"
# We interpret this as printing each term of the calculation.

term1_coeff = (2/3 * lambda1 * math.exp(lambda1 / 2) - 1/3 * math.exp(lambda1 / 2))
term1 = term1_coeff * x2_0

term2 = (2/3) * lambda2 * math.exp(lambda2 / 2)
term3 = (10/3) * math.exp(lambda1 / 2)

# The equation is: Final_Value = term1 - term2 - term3
print(f"Term 1 (coefficient of x2(0) * x2(0)): {term1}")
print(f"Term 2 (value to be subtracted): {term2}")
print(f"Term 3 (value to be subtracted): {term3}")

# Step 5: Calculate the final value
final_value = term1 - term2 - term3
print("\nFinal Equation:")
print(f"{term1} - ({term2}) - ({term3}) = {final_value}")

print("\nResult:")
print(final_value)