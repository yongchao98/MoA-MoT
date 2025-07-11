# Based on the analytical solution, the integral part of the expression evaluates to 5^(1/4) / 12.
# The problem then simplifies to computing (12)^4 * (5^(1/4) / 12)^4.

# Define the numbers that appear in the simplified final equation
base = 12
power = 4
integral_numerator_base = 5
exp_num = 1
exp_den = 4
integral_denominator = 12

# The final result is derived from the simplification:
# (12^4 * 5) / 12^4 = 5
final_result = 5

# Print the full final equation, showing how all the numbers fit together
print(f"({base})^{power} * (({integral_numerator_base}^({exp_num}/{exp_den})) / {integral_denominator})^{power} = {final_result}")