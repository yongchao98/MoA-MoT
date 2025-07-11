def solve_quiver_reflection_questions():
  """
  Solves the three yes/no questions about the quiver with a reflection automorphism.
  The reasoning is provided in the text above.
  """
  
  # Based on the analysis:
  # (a) The action of sigma is expected to mirror the action of g on an arrow at a fixed vertex. Yes.
  # (b) This specific relation is characteristic of underlying algebraic structures like Calabi-Yau algebras. Yes.
  # (c) A counterexample shows this condition is not necessary. The true condition is lambda^2=1. No.

  answer_a = "Yes"
  answer_b = "yes"
  answer_c = "no"

  print(f"({answer_a}); ({answer_b}); ({answer_c}).")

solve_quiver_reflection_questions()