import math

# Step 1: Define the constants lambda1 and lambda2
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# The expression to be calculated is:
# V = ((2/3)*lambda1*e^(lambda1/2) - (1/3)*e^(lambda1/2))*x2(0) - (2/3)*lambda2*e^(lambda2/2) - (10/3)*e^(lambda1/2)
# The value of x2(0) is not given. A standard assumption in controllability problems
# is that the system starts from the origin, so x(0) = 0.
# We will proceed with the assumption that x2(0) = 0.

x2_0 = 0

# With x2(0) = 0, the expression simplifies to:
# V = - (2/3)*lambda2*e^(lambda2/2) - (10/3)*e^(lambda1/2)

# Step 2: Calculate each term of the simplified expression
term1_val = -(2/3) * lambda2 * math.exp(lambda2 / 2)
term2_val = -(10/3) * math.exp(lambda1 / 2)

# Step 3: Calculate the final value
final_value = term1_val + term2_val

# Step 4: Print the calculation as requested
print("Assuming x2(0) = 0, the expression simplifies to: V = - (2/3)*lambda2*e^(lambda2/2) - (10/3)*e^(lambda1/2)")
print("\nLet's plug in the values:")
print(f"lambda1 = {lambda1}")
print(f"lambda2 = {lambda2}")
print("\nThe final equation is:")
print(f"V = -({2/3}) * ({lambda2}) * e^({lambda2 / 2}) - ({10/3}) * e^({lambda1 / 2})")
print(f"\nValue of the first term: {term1_val}")
print(f"Value of the second term: {term2_val}")
print(f"\nFinal calculated value: {final_value}")
