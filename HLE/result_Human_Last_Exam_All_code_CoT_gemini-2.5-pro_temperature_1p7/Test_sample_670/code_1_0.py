# The Tate-style program-expression graph represents the recursive computation
# in the loop. The mu operator is used for this recursion.
# 'a' is the first bound variable for the recursive value of the variable 'a'.
# 'b' is the second bound variable for the induction variable 'i'.
# The operation inside the loop is 'a * b'.

# The final expression is constructed from these parts.
# Note: The expression itself is symbolic and does not contain any numbers
# from the original code, as the initial values are considered inputs to
# the μ-node in the graph, not part of its definition.

final_expression = "μ a b . a * b"

print(final_expression)