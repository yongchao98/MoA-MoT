import math

# The closed form of the infinite product is P = 2^(1/8) * e^(-pi/24).
# This script will display the components of this expression and its numerical value.

# Define the components of the final expression
base1 = 2
exponent1_num = 1
exponent1_den = 8
exponent1 = exponent1_num / exponent1_den

base2_symbol = "e"
exponent2_coeff = -1
exponent2_den = 24

# Construct the symbolic expression as a string
expression_string = f"{base1}^({exponent1_num}/{exponent1_den}) * {base2_symbol}^({exponent2_coeff}*pi/{exponent2_den})"

# Print the final symbolic expression
print("The closed-form expression for the infinite product is:")
print(f"P = {expression_string}")
print("-" * 30)

# Print each number/component in the final equation as requested
print("The components of the final expression are:")
print(f"Term 1 base: {base1}")
print(f"Term 1 exponent: {exponent1_num}/{exponent1_den}")
print(f"Term 2 base: {base2_symbol} (Euler's number)")
print(f"Term 2 exponent: ({exponent2_coeff} * pi) / {exponent2_den}")
print("-" * 30)

# Calculate the numerical value for approximation
pi_val = math.pi
numerical_value = (base1**exponent1) * math.exp((exponent2_coeff * pi_val) / exponent2_den)

# Print the numerical value
print(f"The numerical value of the expression is approximately:")
print(numerical_value)