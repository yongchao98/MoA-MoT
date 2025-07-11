def solve_group_theory_question():
  """
  This function provides the solution to the given group theory problem.
  The reasoning is provided in the text above.
  """

  # Part (a): Existence and uniqueness of the minimal group hat(G).
  # The described group hat(G) is the p-localization of G. Standard theorems
  # in group theory guarantee its existence and uniqueness up to a G-isomorphism.
  answer_a = "Yes"

  # Part (b): Maximum possible derived length of hat(G).
  # As derived in the thinking steps, the derived length of hat(G) is bounded
  # by n. This bound can be achieved by specific constructions of the group G.
  # The maximum possible derived length is therefore n.
  # The problem asks for a single numerical value, but since n is a parameter,
  # the answer is an expression in terms of n.
  answer_b_expression = "n"

  # The question asks for the output in a specific format.
  # This code prints the answer string.
  # Final answer format required: (a) [Yes/No]; (b) [Expression].
  final_answer = f"(a) {answer_a}; (b) {answer_b_expression}"
  print(final_answer)

solve_group_theory_question()