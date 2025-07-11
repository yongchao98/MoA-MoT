def solve_complexity_questions():
  """
  Prints the answers to the complexity theory questions based on the analysis.
  """
  
  # (a) If `IS_NOT_FREE` is NP-hard, does that imply `IS_FREE` is also NP-hard?
  # Analysis: No. It implies `IS_FREE` is coNP-hard.
  answer_a = "No"

  # (b) If `IS_NOT_FREE` is NP-complete, does that imply `IS_FREE` is also NP-complete?
  # Analysis: No. It implies `IS_FREE` is coNP-complete.
  answer_b = "No"

  # (c) If `IS_FREE` is in NP, and `IS_NOT_FREE` is NP-hard, does that imply `IS_FREE` is NP-complete?
  # Analysis: Yes. The premises together imply NP = coNP, which in turn makes `IS_FREE` NP-hard.
  # Since it is also in NP (by premise), it is NP-complete.
  answer_c = "Yes"

  print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_complexity_questions()