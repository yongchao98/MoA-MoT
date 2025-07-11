import math

# The final result for y(0) is derived from the expression (3/4) * (9/4)^(3/5).
# Here, we define the individual numbers that make up this final equation.
numerator_outer = 3
denominator_outer = 4
base_numerator = 9
base_denominator = 4
exp_numerator = 3
exp_denominator = 5

# We print the final equation in a readable format, showing each number.
print("The final equation for the deflection y(0) is:")
print(f"y(0) = ({numerator_outer}/{denominator_outer}) * ({base_numerator}/{base_denominator})^({exp_numerator}/{exp_denominator})")
print("-" * 30)

# Now, we perform the calculation step-by-step.
# Calculate the value of the base.
base_val = base_numerator / base_denominator

# Calculate the value of the exponent.
exponent_val = exp_numerator / exp_denominator

# Calculate the base raised to the power of the exponent.
powered_base_val = base_val ** exponent_val

# Calculate the final result for y(0).
y_0 = (numerator_outer / denominator_outer) * powered_base_val

# Print the final calculated deflection.
print(f"The calculated value for y(0) is: {y_0}")