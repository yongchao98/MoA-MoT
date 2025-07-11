def solve_group_theory_question():
  """
  This function provides the solution to the group theory problem.
  
  The solution is based on the following reasoning:
  
  Part (a): The group hat(G) is the p-completion of G. Standard theory on group completions
  confirms that such a minimal group exists and is unique up to isomorphism.
  
  Part (b): The derived length of hat(G) is the same as the derived length of G.
  The given subnormal series is G = G_1 > G_2 > ... > G_n > G_{n+1}={1}, where each
  factor G_i/G_{i+1} is abelian. This implies G_i' is a subgroup of G_{i+1}.
  An induction argument shows that the k-th derived subgroup G^(k) is a subgroup of G_{k+1}.
  For k=n, this means G^(n) is a subgroup of G_{n+1}={1}, so dl(G) <= n.
  This maximum is attainable, for example by the iterated wreath product of n copies of
  a cyclic p-group, which has derived length n and satisfies the series condition.
  Thus, the maximum possible derived length is n.
  """
  
  answer_a = "Yes"
  
  # The expression for the maximum derived length is the length of the series, n.
  # The question asks for an expression, which in this case is a variable from the prompt.
  answer_b_expression = "n"
  
  # Final answer in the specified format
  final_answer = f"(a) {answer_a}; (b) {answer_b_expression}"
  
  print(final_answer)

solve_group_theory_question()