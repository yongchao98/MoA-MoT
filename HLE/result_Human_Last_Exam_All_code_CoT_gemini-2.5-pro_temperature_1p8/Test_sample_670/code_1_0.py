import unicodedata

# This program generates the Tate-style program-expression for the variable 'a'
# as defined within the provided loop.

# The core of the problem is to represent the recursive definition of the loop's state.
# The state is determined by the variables 'a' and 'i'.
# In the Tate framework, this recursion is captured by the μ (mu) operator.

# 1. Define the bound variables. The prompt specifies using 'a' for the first
#    (representing the program variable 'a') and 'b' for the second
#    (representing the program variable 'i').
bound_vars = "(a, b)"

# 2. Define the update expressions for the state variables in terms of the bound variables.
#    - The new value of 'a' is the old 'a' times the old 'i', which is 'a * b'.
#    - The new value of 'i' is the old 'i' plus one, which is 'b + 1'.
#    The number '1' is explicitly part of the final expression.
update_a_expr = "a * b"
constant_number = 1
update_b_expr = f"b + {constant_number}"

# 3. Assemble the full μ-expression.
#    The structure is: μ (bound_vars) . (update_expr_for_a, update_expr_for_b)
mu_char = unicodedata.lookup("GREEK SMALL LETTER MU")
final_expression = f"{mu_char} {bound_vars} . ({update_a_expr}, {update_b_expr})"

# 4. Print the final expression to the console.
print(final_expression)