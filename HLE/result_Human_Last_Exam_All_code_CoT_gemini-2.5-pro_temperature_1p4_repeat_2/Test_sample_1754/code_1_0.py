import math

# Define the given lambda values
lambda1 = (1 + math.sqrt(17)) / 2
lambda2 = (1 - math.sqrt(17)) / 2

# From the analysis, we assume the system is in an equilibrium state
# that satisfies the boundary condition. This leads to x_2(0) = 0.
x2_0 = 0

# The expression to calculate is:
# (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2(0) 
# - 2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)

# Substitute x2(0) = 0 into the expression
term1_coeff_val = (2/3 * lambda1 * math.exp(lambda1/2) - 1/3 * math.exp(lambda1/2))
term1_val = term1_coeff_val * x2_0

term2_val = (2/3) * lambda2 * math.exp(lambda2/2)
term3_val = (10/3) * math.exp(lambda1/2)

# Calculate the final value
final_value = term1_val - term2_val - term3_val

# Print the equation with all the numbers
print(f"Given lambda1 = {lambda1}")
print(f"Given lambda2 = {lambda2}")
print(f"From the problem analysis, we deduce that x2(0) = {x2_0}")
print("\nCalculating the expression:")
print(f"({2/3}*({lambda1})*exp({lambda1}/2) - {1/3}*exp({lambda1}/2))*({x2_0}) - {2/3}*({lambda2})*exp({lambda2}/2) - {10/3}*exp({lambda1}/2)")
print(f"= ({term1_coeff_val})*({x2_0}) - {term2_val} - {term3_val}")
print(f"= {term1_val} - {term2_val} - {term3_val}")
print(f"= {final_value}")

print("\nFinal Answer:")
print(f"{final_value}")