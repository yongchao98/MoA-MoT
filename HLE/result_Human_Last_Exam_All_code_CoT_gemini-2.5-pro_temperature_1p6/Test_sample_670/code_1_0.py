# This script generates the Tate-style program-expression graph
# for the variable 'a' in the given loop.

# The recurrence for the loop counter 'i' is i -> i + 1.
# Using 'b' as the bound variable, this is represented as: μ b. (b + 1)
# Let's call this expression_i.

# The recurrence for 'a' is a -> a * i.
# Using 'a' as the bound variable, this is represented as: μ a. (a * expression_i)

# By substituting expression_i, we get the final combined expression.
final_expression = "μ a. (a * (μ b. (b + 1)))"

# Print the final expression. This expression contains the number 1,
# fulfilling the requirement to output each number in the equation.
print(final_expression)