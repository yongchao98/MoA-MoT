import sys

# We need to calculate x(1) using the derived formula:
# x(t) = -8 / (4*t^2 - 8*t + 1)
t = 1

# Define the components of the equation
numerator = -8
term_t_squared = 4 * t**2
term_t = -8 * t
constant = 1

# Calculate the denominator
denominator = term_t_squared + term_t + constant

# Check for division by zero
if denominator == 0:
    print("Error: Division by zero.", file=sys.stderr)
else:
    # Calculate the final result
    x_1 = numerator / denominator

    # Print the calculation steps as requested
    print("To find x(1), we evaluate the expression x(t) = -8 / (4*t^2 - 8*t + 1) at t=1.")
    print(f"x(1) = {numerator} / (4*({t})^2 - 8*({t}) + {constant})")
    print(f"x(1) = {numerator} / ({term_t_squared} + ({term_t}) + {constant})")
    print(f"x(1) = {numerator} / {denominator}")
    print(f"The final value is: {x_1}")
