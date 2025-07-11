def solve_graph_theory_problem():
  """
  This function formats and prints the final answer based on the analysis.
  """
  
  # Part (a) Answer: True
  # As shown in the analysis, the function is bounded both above and below by linear functions of n.
  answer_a = "True"
  
  # Part (b) Answer: True
  # The number of edges is bounded by a constant that depends only on s and t, not n.
  answer_b = "True"
  
  # Part (c) Answer: (s-1)(2s + 2t - 5)
  # This expression for the upper bound was derived from the vertex cover argument.
  # The format below respects the instruction to "output each number".
  expression_c = "(s - 1) * (2 * s + 2 * t - 5)"
  
  # The final formatted string
  final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
  
  print(final_answer)

solve_graph_theory_problem()