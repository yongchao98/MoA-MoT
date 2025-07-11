def solve_extremal_problem():
  """
  This function provides the solution to the Turan-type extremal problem.
  """
  
  # Part (a)
  # If G is not a union of K_2's, ex(n; G, K_{1,t}-ind) = Theta(n).
  # This is True as explained in the analysis.
  answer_a = "True"
  
  # Part (b)
  # If G ~ sK_2, ex(n; sK_2, K_{1,t}-ind) is independent of n.
  # This is True as the number of edges is bounded by a constant depending on s and t.
  answer_b = "True"
  
  # Part (c)
  # The upper bound for ex(n; sK_2, K_{1,t}-ind) in terms of s and t.
  # The expression is (s-1)(2s + 2t - 5).
  s = "s"
  t = "t"
  # The numbers in the equation are 1, 2, 2, 5.
  answer_c = f"({s}-1)*(2*{s} + 2*{t} - 5)"
  
  print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_extremal_problem()