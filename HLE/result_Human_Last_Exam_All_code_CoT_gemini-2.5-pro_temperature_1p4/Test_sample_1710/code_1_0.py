def solve_go_problem():
  """
  This function provides the solution to the Go problem.
  Based on the analysis, there are two moves that successfully
  capture the two white stones and ensure the black group's survival.
  """
  # The optimal moves are to play at (2, 5) or (3, 5).
  # These moves create a "net" (geta) that makes the capture of the
  # two white stones at (2, 6) and (3, 6) inevitable.
  # The answers are ordered lexicographically as requested.
  answer1 = (2, 5)
  answer2 = (3, 5)

  print(f"({answer1[0]}, {answer1[1]}), ({answer2[0]}, {answer2[1]})")

solve_go_problem()