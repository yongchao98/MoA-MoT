def solve_turán_problem():
  """
  This function provides the solution to the Turán-type extremal problem.
  It constructs and prints the final answer string based on the reasoning above.
  """

  # Part (a): True. The extremal function is linear in n.
  answer_a = "True"

  # Part (b): True. The extremal function is bounded by a constant, independent of n for large n.
  answer_b = "True"

  # Part (c): The derived upper bound expression.
  # The prompt asks to output each number in the final equation.
  # The expression is (s-1)(2s+2t-5). The numbers are 1, 2, 2, 5.
  expression_c = "(s - 1)*(2*s + 2*t - 5)"

  # Combine the answers into the specified format.
  final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"

  print(final_answer)

solve_turán_problem()