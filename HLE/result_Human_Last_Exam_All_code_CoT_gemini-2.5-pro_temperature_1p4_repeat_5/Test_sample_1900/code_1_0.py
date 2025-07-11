def get_complexity_answer():
  """
  This function provides the big-O notation for the computational complexity
  of the two questions asked.

  Question A: The complexity of deciding if a line exists through all n red balls.
  This is the Hamiltonian Path decision problem, which is NP-complete. A standard
  exact algorithm runs in O(n^2 * 2^n).

  Question B: The complexity of finding such a line, given it exists.
  This is the search version of the Hamiltonian Path problem. Its complexity is
  asymptotically the same as the decision problem.
  """
  complexity_A = "O(n^2 * 2^n)"
  complexity_B = "O(n^2 * 2^n)"

  # The problem asks for the two answers to be separated by a semicolon.
  print(f"{complexity_A}; {complexity_B}")

get_complexity_answer()