# The value of 'a' in the loop is determined by a recurrence relation involving
# both 'a' and the loop counter 'i'. We can express these recurrences using
# the μ (mu) operator from Tate's program-expression graphs.
# The Φ (phi) function represents the merging of values from different control
# paths (i.e., the initial entry into the loop and the subsequent iterations).

# 1. Expression for the loop counter 'i':
#    - It starts at 1.
#    - In each iteration, it's incremented by 1.
#    - Using 'b' as the bound variable, the expression is: μ b . Φ(1, b + 1)
#      where 1 is the initial value and b + 1 is the recurrence.

# 2. Expression for the variable 'a':
#    - It is initialized to 1 before the loop.
#    - In each iteration, it's updated by the rule: a = a * i.
#    - Using 'a' as the bound variable, its expression depends on the expression for 'i'.
#    - The expression is: μ a . Φ(1, a * [expression_for_i])

# 3. Combining these gives the final nested expression.

# The following Python code constructs and prints this expression.

a_initial_val = 1
i_initial_val = 1
i_increment = 1

# Using f-string formatting to embed the numbers directly into the final expression string.
# The inner expression represents 'i'.
inner_expression = f"μ b . Φ({i_initial_val}, b + {i_increment})"

# The outer expression represents 'a' and contains the inner expression.
final_expression = f"μ a . Φ({a_initial_val}, a * ({inner_expression}))"

print(final_expression)