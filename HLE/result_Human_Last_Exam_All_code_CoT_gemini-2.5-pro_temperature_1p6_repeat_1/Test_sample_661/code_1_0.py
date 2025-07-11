import sys

# The Knuth-Bendix completion algorithm is a complex procedure.
# The following code prints the result derived from applying the algorithm
# to the given term-rewriting system, as explained in the steps above.

# The added rules are:
# 1. f(g(x), g(y)) -> g(x)
# 2. g(g(g(x))) -> g(x)
# 3. h(x) -> g(x)
#
# These rules are ordered increasingly by their left-hand side
# using the lexicographic path ordering induced by f < g < h.

def solve():
  """
  Lists the rules added by the Knuth-Bendix completion algorithm for the given system.
  """
  # The rules are ordered increasingly by LHS using the LPO f<g<h.
  # LHS ordering: f(g(x), g(y)) < g(g(g(x))) < h(x)
  rule1 = "f(g(x), g(y)) -> g(x)"
  rule2 = "g(g(g(x))) -> g(x)"
  rule3 = "h(x) -> g(x)"

  # The output is a comma-separated list of the new rules.
  final_answer = f"{rule1}, {rule2}, {rule3}"
  print(final_answer)

solve()