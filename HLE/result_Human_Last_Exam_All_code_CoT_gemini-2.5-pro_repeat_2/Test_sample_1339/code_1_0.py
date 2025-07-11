def solve_group_theory_question():
  """
  This function provides the answer to the posed group theory question.
  """

  # Part (a): Existence and uniqueness of the minimal group.
  # Based on the theory of completions of solvable groups, the answer is Yes.
  answer_a = "Yes"

  # Part (b): Maximum possible derived length.
  # As explained in the reasoning, the derived length of hat(G) is the same as G.
  # The subnormal series of length n implies dl(G) <= n.
  # This bound is achievable.
  # The maximum possible derived length is therefore n.
  answer_b_expression = "n"

  # Print the final answer in the requested format.
  print(f"(a) [{answer_a}]; (b) [{answer_b_expression}]")

solve_group_theory_question()