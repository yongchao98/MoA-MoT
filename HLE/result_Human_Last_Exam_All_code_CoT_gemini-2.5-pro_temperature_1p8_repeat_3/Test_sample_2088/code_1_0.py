# The original problem is to compute (12)^4 * (integral)^4
# After simplification and evaluation, the value of the integral is found to be 5**(1/4) / 12.

# Define the components of the expression
base = 12
exponent = 4
integral_value_numerator_base = 5
integral_value_numerator_exponent = 0.25
integral_value_denominator = 12

# The integral value as a symbolic representation
integral_val_str = f"{integral_value_numerator_base}**{integral_value_numerator_exponent} / {integral_value_denominator}"

print(f"The expression to compute is ({base})^{exponent} * (I_1 - I_2)^{exponent}")
print(f"The value of the integral part, I_1 - I_2, is {integral_val_str}")

# The full expression becomes:
# (12)^4 * ( (5**(1/4))/12 )^4
# This simplifies to 12^4 * (5 / 12^4) which is 5.
final_result = (base ** exponent) * ((integral_value_numerator_base ** integral_value_numerator_exponent / integral_value_denominator) ** exponent)

print("\nThe calculation is:")
print(f"({base})^{exponent} * ({integral_val_str})^{exponent}")
print(f"= ({base**exponent}) * (({integral_value_numerator_base}^{integral_value_numerator_exponent*exponent}) / ({integral_value_denominator}^{exponent}))")
print(f"= ({base**exponent}) * ({int(integral_value_numerator_base**(integral_value_numerator_exponent*exponent))} / {integral_value_denominator**exponent})")
print(f"= {int(final_result)}")
