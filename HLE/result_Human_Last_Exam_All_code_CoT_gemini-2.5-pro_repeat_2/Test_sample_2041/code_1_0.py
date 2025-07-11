# The number of distinct functions is determined by the number of boolean functions
# that can be constructed from the available shallow atomic expressions.

# We found there are 4 independent atomic expressions that can be formed:
# p(q_x), p(q_not_x), p(q_T), and p(q_F).
# These can be thought of as 4 boolean variables.
num_vars = 4

# The total number of distinct functions is the number of boolean functions of `num_vars` variables.
# For n variables, this number is 2^(2^n).
base1 = 2
base2 = 2
exponent = num_vars

# Calculate the result
result = base1 ** (base2 ** exponent)

# Print the equation and the final answer
print(f"{base1} ** ({base2} ** {exponent}) = {result}")
