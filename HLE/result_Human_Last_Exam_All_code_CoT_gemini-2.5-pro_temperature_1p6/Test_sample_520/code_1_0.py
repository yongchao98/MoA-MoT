def find_profession():
  """
  This function provides the answer to the historical question about Andre Vernault.
  """
  profession = "lion-tamer"
  explanation = (
      "During World War II, counterintelligence in the UK falsely believed Belgian refugee Andre Vernault to be a German spy "
      "because he concealed his true profession from interrogators. He was embarrassed to admit that he was a {}."
  ).format(profession)
  print(explanation)

find_profession()