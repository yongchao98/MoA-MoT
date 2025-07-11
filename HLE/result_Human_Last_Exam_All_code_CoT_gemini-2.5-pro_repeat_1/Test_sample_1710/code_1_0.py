def solve_go_problem():
  """
  This function determines the correct moves for Black to survive and
  capture the two white stones on the given Go board.

  The analysis identifies two vital points that achieve both objectives.
  The points are (2, 5) and (3, 5).
  """
  
  # The two correct moves for Black, in lexicographical order.
  # Row is counted from top to bottom (1-8), column from left to right (1-9).
  answer_1_row = 2
  answer_1_col = 5
  answer_2_row = 3
  answer_2_col = 5

  # Print the result in the specified format.
  print(f"({answer_1_row}, {answer_1_col}), ({answer_2_row}, {answer_2_col})")

solve_go_problem()