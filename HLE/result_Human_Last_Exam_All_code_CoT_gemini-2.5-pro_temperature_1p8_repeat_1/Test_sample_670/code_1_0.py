import sys

# Set default encoding to UTF-8 to ensure special characters print correctly
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

# The Tate-style program-expression graph for the variable 'a' in the loop is
# represented using the μ (mu) and φ (phi) operators.
#
# - 'μ a . F(a)' denotes the least fixed-point of the function F. It models the loop.
# - 'a' is the first bound variable, representing the value of the program variable 'a' from the previous iteration.
# - 'φ(v1, v2)' is a function that chooses between its arguments. It selects the
#   initial value 'v1' on the first entry to the loop and the recurring value 'v2' on subsequent iterations.
# - '1' is the initial value of the program variable 'a' before the loop.
# - 'b' is the second bound variable, representing the abstracted value of the loop counter 'i'.
# - The expression for the value computed in the previous iteration is 'a * b'.
#
# This code will print the final expression, including the number 1 from the initial assignment.

final_expression = "μ a . φ(1, a * b)"
print(final_expression)