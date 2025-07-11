# Define the components of the expression from the simplified form
base = 12
power = 4
integral_value_numerator = 5

# The expression to compute is (base**power) * I**power
# I = (integral_value_numerator**(1/4)) / base
# I**power = integral_value_numerator / (base**power)
I_pow_4 = integral_value_numerator / (base**power)
result = (base**power) * I_pow_4

# Print the final equation with all numbers
print("Let I be the result of the combined integrals after substitution.")
print(f"The value of the integral I is ({integral_value_numerator}^(1/4)) / {base}.")
print("The expression to compute is (12)^4 * I^4.")
print(f"I^4 = (({integral_value_numerator}^(1/4)) / {base})^4 = {integral_value_numerator} / ({base}^{power}).")
print(f"The final computation is: ({base})^{power} * ({integral_value_numerator} / {base}^{power}) = {int(result)}")