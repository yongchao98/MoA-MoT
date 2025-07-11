def solve_dessin_question():
  """
  Solves the three-part question about dessin d'enfants and prints the answer.
  """
  answer_a = "Yes"
  answer_b = "HA, TW, AS"
  answer_c = "False"

  final_answer = f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]."
  print(final_answer)
  # Adding a final line to be captured by the special format string.
  # The problem does not have a numerical answer, so we format the text answer itself.
  print("<<<" + f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]." + ">>>")

solve_dessin_question()