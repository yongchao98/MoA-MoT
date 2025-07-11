def solve():
  """
  This function provides the answer for the computational complexity questions.
  """
  # Answer for Question A: Complexity of deciding whether a line through all red balls exists.
  # The properties of the graph (connected and locally connected) guarantee that a Hamiltonian path always exists.
  # Therefore, the decision problem is trivial and takes constant time.
  answer_A = "O(1)"

  # Answer for Question B: Complexity of finding a line that goes through all the red balls.
  # This is the search problem for a Hamiltonian path in a connected, locally connected graph.
  # While the general problem is NP-hard, for this specific class of graphs, polynomial-time algorithms exist.
  # Known algorithms have complexities such as O(n^5), O(n^4), or O(n^3). We provide O(n^3) as a
  # representative complexity for a known polynomial-time algorithm.
  answer_B = "O(n^3)"

  # The final output should be in the format O(ans_A); O(ans_B)
  print(f"{answer_A}; {answer_B}")

solve()