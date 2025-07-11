def solve_extremal_problem():
  """
  Solves the given graph theory problem and prints the answer.
  """

  # Part (a)
  # The statement is false. A counterexample is any non-bipartite graph G
  # that is not a matching, such as G = C_3. For such a G, the extremal
  # function is known to be super-linear, not Theta(n).
  answer_a = "False"

  # Part (b)
  # The statement is true. The number of edges in a graph that is sK2-free
  # and induced K_{1,t}-free is bounded by a constant that depends on s and t,
  # but not on n. A non-decreasing integer function bounded above must be
  # constant for large enough n.
  answer_b = "True"

  # Part (c)
  # An upper bound can be derived by partitioning the graph H into V(M)
  # (vertices of a max matching) and I (an independent set).
  # Edges <= edges_in_V(M) + edges_between_V(M)_and_I
  # Edges <= C(2(s-1), 2) + 2(s-1)*(t-1)
  # Simplifying this expression gives (s-1)*(2s + 2t - 5).
  # The numbers in the equation are 1, 2, 2, 5.
  expression_c = "(s - 1) * (2 * s + 2 * t - 5)"

  # Format the final output
  final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
  print(final_answer)

solve_extremal_problem()