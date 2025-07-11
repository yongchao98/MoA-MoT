# The original expression to compute is (12)^4 * (Integral_1 - Integral_2)^4.
# Based on the step-by-step analytical derivation, the value of the difference of the integrals is (5^(1/4))/12.
# Here we calculate the final value and print the equation.

# Define the numbers from the expression and the result of the integral
val_12 = 12
power_4 = 4
integral_numerator_base = 5

# Calculate the final result
# The expression becomes (12^4) * ( (5^(1/4) / 12)^4 )
# This simplifies to (12^4) * ( 5 / 12^4 ), which is 5.
result = (val_12 ** power_4) * (((integral_numerator_base**(1/power_4))/val_12)**power_4)

# Print the final equation with each number shown
# The string formatting will display the equation used for the final calculation based on the derived integral value.
# f-string formatting for 1/4 results in 0.25, so we represent it as (1/4) for clarity.
print(f"({val_12})^{power_4} * (({integral_numerator_base}^(1/{power_4}))/{val_12})^{power_4} = {int(round(result))}")