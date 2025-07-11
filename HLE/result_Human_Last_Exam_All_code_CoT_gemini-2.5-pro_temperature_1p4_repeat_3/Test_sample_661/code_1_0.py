def solve_knuth_bendix():
  """
  This function executes the plan to find the new rules added by the Knuth-Bendix completion
  and prints them in the specified order.
  """

  # These are the rules identified as being added during the completion process.
  rule1_lhs = "f(g(x), g(y))"
  rule1_rhs = "g(x)"

  rule2_lhs = "g(g(g(x)))"
  rule2_rhs = "g(x)"

  rule3_lhs = "h(x)"
  rule3_rhs = "g(x)"

  # The rules are ordered increasing by their Left-Hand-Side using LPO with f < g < h.
  # Order: f(g(x), g(y)) < g(g(g(x))) < h(x)
  ordered_rules = [
    f"{rule1_lhs} -> {rule1_rhs}",
    f"{rule2_lhs} -> {rule2_rhs}",
    f"{rule3_lhs} -> {rule3_rhs}"
  ]

  # The output is a comma-separated string of the rules.
  print(", ".join(ordered_rules))

solve_knuth_bendix()