import math

# Define the components of the expression from the problem
base = 12
exponent = 4

# Based on our analytical derivation, the value of the integral is (5^(1/4))/12
# Let's define the components for the value of the integral
integral_numerator_base = 5
integral_denominator = 12

# Calculate the value of the integral
integral_value = integral_numerator_base**(1/exponent) / integral_denominator

# Calculate the final result of the expression (base^exponent) * (integral_value^exponent)
result = base**exponent * integral_value**exponent

# To meet the request of outputting each number in the final equation,
# we format the print statement to show the structure of the calculation.
# The final equation is (12)^4 * ( (5^(1/4))/12 )^4 = 5.
print("The value of the expression is computed as follows:")
print(f"({base})^({exponent}) * ({integral_numerator_base}^(1/{exponent}) / {integral_denominator})^({exponent}) = {int(round(result))}")