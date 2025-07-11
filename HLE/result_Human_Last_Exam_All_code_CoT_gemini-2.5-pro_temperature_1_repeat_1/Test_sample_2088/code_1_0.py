# The original expression is (12)^4 * (Integral)^4
# As derived in the explanation, the value of the integral is (5^(1/4)) / 12.

# Define the numbers from the problem and our derivation
base = 12
power = 4
integral_numerator_base = 5
integral_denominator = 12

# Calculate the components of the final expression
term1 = base ** power
integral_val_numerator = integral_numerator_base ** (1/4)
integral_val_denominator = integral_denominator
integral_val_raised_to_power_numerator = (integral_val_numerator) ** power
integral_val_raised_to_power_denominator = (integral_val_denominator) ** power

# Calculate the final result
final_result = term1 * (integral_val_raised_to_power_numerator / integral_val_raised_to_power_denominator)

# Print the equation with all the numbers
print("The calculation is:")
print(f"({base})^{power} * (({integral_numerator_base}^(1/4)) / {integral_denominator})^{power} = {term1} * ({int(integral_val_raised_to_power_numerator)} / {integral_val_raised_to_power_denominator}) = {int(final_result)}")