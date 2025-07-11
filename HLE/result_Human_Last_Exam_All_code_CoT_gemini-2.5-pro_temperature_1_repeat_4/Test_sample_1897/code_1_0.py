# Define the numbers in the final expression
base = 2
exponent_numerator = 1
exponent_denominator = 8
pi_denominator = 24

# Construct and print the final closed expression
# The expression is of the form: base^(exponent_numerator/exponent_denominator) * e^(-pi/pi_denominator)
final_expression = f"{base}^({exponent_numerator}/{exponent_denominator}) * e^(-pi/{pi_denominator})"

print("The closed form of the infinite product is:")
print(final_expression)