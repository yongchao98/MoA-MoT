def solve_closepact():
  """
  Determines which of the given sets are necessarily "closepact".

  Based on the analysis, a subset of R or C is closepact if and only if it is compact.
  We identify the compact sets from the list.

  C. A finite set is compact.
  G. A convergent sequence union its limit is compact.
  H. A convergent sequence union its limit is compact.
  J. A closed interval is compact.
  M. The Cantor set is compact.

  The final answer is the concatenation of these letters.
  """
  answer = "CGHJM"
  print(answer)

solve_closepact()