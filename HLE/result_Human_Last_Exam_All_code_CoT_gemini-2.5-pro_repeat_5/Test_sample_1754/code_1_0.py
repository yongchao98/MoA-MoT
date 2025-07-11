import math

# Step 1: Define the given constants lambda_1 and lambda_2
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# Step 2: The problem is ambiguous as x2(0) is not given.
# We assume the problem starts from the origin, a common case in controllability, so x(0) = 0, which means x2(0) = 0.
# The expression to calculate V simplifies to:
# V = - (2/3) * lambda2 * exp(lambda2 / 2) - (10/3) * exp(lambda1 / 2)
x2_0 = 0

# Step 3: Calculate each part of the expression.
# The first major term in the original expression evaluates to zero with our assumption.
term1_full_expr = ((2/3) * lambda1 * math.exp(lambda1 / 2) - (1/3) * math.exp(lambda1 / 2)) * x2_0

# The remaining terms are term2 and term3
term2 = (2/3) * lambda2 * math.exp(lambda2 / 2)
term3 = (10/3) * math.exp(lambda1 / 2)

# The final value is V = term1_full_expr - term2 - term3
final_value = term1_full_expr - term2 - term3

# Step 4: Print the equation with the calculated numbers.
# The problem asks to find the value of:
# (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2(0) - 2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)
# With x2(0) = 0, this is: 0 - term2 - term3

print("Based on the assumption that the process starts from the origin, x2(0) = 0.")
print("The expression to calculate is:")
print(f"({(2/3)*lambda1*math.exp(lambda1/2) - (1/3)*math.exp(lambda1/2):.4f}) * {x2_0} - ({term2:.4f}) - ({term3:.4f})")
print("\nThis simplifies to:")
print(f"0 - {term2:.4f} - {term3:.4f} = {final_value:.4f}")
print(f"\nThe final calculated value is: {final_value:.4f}")