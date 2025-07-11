# The initial value problem is x'(t) = (t-1)x^2(t), with x(0) = -8.
# By separating variables and integrating, we get:
# -1/x = t^2/2 - t + C
# Using the initial condition x(0) = -8, we find C = 1/8.
# So, the particular solution is -1/x = t^2/2 - t + 1/8.
# Solving for x(t) gives: x(t) = -8 / (4*t^2 - 8*t + 1).

# We need to find x(1). Let's define the variable and the numbers.
t = 1
numerator = -8
coeff_t_squared = 4
coeff_t = -8
constant = 1

# Calculate the denominator of the expression for x(1)
denominator = coeff_t_squared * (t**2) + coeff_t * t + constant

# Calculate the final result
result = numerator / denominator

# Print out the calculation steps as requested
print(f"The particular solution is x(t) = -8 / (4*t^2 - 8*t + 1).")
print(f"To find x(1), we substitute t = 1 into this equation:")
print(f"x(1) = {numerator} / ({coeff_t_squared} * ({t})**2 - {abs(coeff_t)} * ({t}) + {constant})")
print(f"x(1) = {numerator} / ({coeff_t_squared * t**2} - {abs(coeff_t) * t} + {constant})")
print(f"x(1) = {numerator} / ({denominator})")
print(f"The final value is:")
print(result)