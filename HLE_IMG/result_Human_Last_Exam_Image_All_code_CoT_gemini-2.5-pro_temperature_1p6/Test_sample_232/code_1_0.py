def solve_disease_riddle():
  """
  Identifies the atrophied organ part and the disorder from the provided image context.
  """
  atrophied_organ_part = "caudate nucleus"
  disorder_name = "Huntington's disease"

  # The requested format is "organ, disorder"
  final_answer = f"{atrophied_organ_part}, {disorder_name}"
  print(final_answer)

solve_disease_riddle()