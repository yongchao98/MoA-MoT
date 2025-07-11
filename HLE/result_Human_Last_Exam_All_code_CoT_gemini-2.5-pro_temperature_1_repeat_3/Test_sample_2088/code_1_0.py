# Based on the analytical derivation, the value of the integral difference is (5^(1/4))/12.
# We will use this value to compute the final expression.

# The numbers from the original expression
term1_base = 12
term_power = 4

# The numbers from the derived value of the integral
integral_numerator_base = 5
integral_numerator_power_num = 1
integral_numerator_power_den = 4
integral_denominator = 12

# Calculate the final result
# E = (term1_base^term_power) * ((integral_numerator_base^(1/4) / integral_denominator)^term_power)
result = (term1_base ** term_power) * (((integral_numerator_base ** (integral_numerator_power_num / integral_numerator_power_den)) / integral_denominator) ** term_power)

# Print the final equation with all the numbers
print(f"({term1_base})^{term_power} * (({integral_numerator_base}^({integral_numerator_power_num}/{integral_numerator_power_den})) / {integral_denominator})^{term_power} = {result}")