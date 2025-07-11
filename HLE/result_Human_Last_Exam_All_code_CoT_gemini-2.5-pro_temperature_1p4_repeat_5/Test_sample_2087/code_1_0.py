# The task is to provide the mathematical expression for the limiting CDF.
# This script constructs and prints the formula as a string.

# Define the names of the mathematical terms as strings
term_x = "x"
term_cdf = "F_{X_i}(x)"
term_integral = "I_{X_i}(x)"
term_mean = "mu_{X_i}"

# Build the numerator of the expression
numerator_str = f"({term_x} * {term_cdf} - {term_integral})"

# Build the full expression string
# The final equation is the numerator divided by the mean.
final_expression = f"{numerator_str} / {term_mean}"

# Print the final result to the console
print("The expression for lim_{t->inf} F_{X(t)}(x) is:")
print(final_expression)
