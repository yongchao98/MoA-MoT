def solve_and_print_answer():
  """
  This function determines and prints the solution to the given group theory problem.
  The variable 'n' is treated as a symbolic parameter from the problem statement.
  """

  # Part (a): Existence and uniqueness of hat(G)
  # The group hat(G) as described is the p-localization of G. Its existence and
  # uniqueness (up to isomorphism) are standard results in group theory,
  # established via a universal mapping property.
  answer_a = "Yes"

  # Part (b): Maximum possible derived length of hat(G)
  # The derived length of G is at most n. For the class of groups specified,
  # the p-localization functor preserves the derived length in the maximal case.
  # An upper bound dl(hat(G)) <= n can be established. This maximum is achievable,
  # for instance, by the n-th iterated wreath product of Z with itself, whose
  # derived length is n and whose p-localization also has derived length n.
  # Thus, the maximum possible derived length is n.
  answer_b = "n"

  # The problem asks for the answer in the format (a) [Yes/No]; (b) [Expression].
  # The expression for part (b) is 'n'. There is no equation with numbers to output
  # as the answer is symbolic.
  print(f"(a) {answer_a}; (b) {answer_b}")

solve_and_print_answer()