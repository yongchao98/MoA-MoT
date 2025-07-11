import sys

def solve():
  """
  Solves the theoretical group theory problem.
  """

  # Part (a): Analysis on the existence and uniqueness of the minimal group hat(G).
  # The properties described for hat(G) are those of the p-localization of G.
  # The p-localization construction is universal and yields a unique group up to
  # a G-isomorphism.
  answer_a = "Yes"

  # Part (b): Analysis on the maximum possible derived length of hat(G).
  # The group G has a subnormal series of length n with abelian factors,
  # which implies G is solvable with derived length dl(G) <= n.
  # The p-localization process preserves the derived length for solvable groups,
  # so dl(hat(G)) = dl(G).
  # The maximum value n is attainable, for example, by constructing G as an
  # iterated wreath product of n copies of the integers.
  # The variable n is given symbolically in the problem.
  answer_b = "n"

  # Print the final answer in the required format.
  # The problem states "output each number in the final equation!".
  # Our derived equation is dl(hat(G))_max = n.
  # We will output n as the expression for part (b).
  print(f"(a) {answer_a}; (b) {answer_b}")

solve()