def find_foreign_symbol():
  """
  Analyzes the image and provides an answer to the user's question.
  """
  is_foreign_symbol_present = "True"
  symbol_identity = "modern alphabetic graffiti appearing to be the letters 'no'"
  explanation = (
      f"{is_foreign_symbol_present}. "
      f"In the upper right portion of the image, there is {symbol_identity}. "
      "This writing system was brought to the Americas by Europeans and is unrelated to the pre-Columbian cultures that created the rock art. "
      "Therefore, it is a symbol that is not from the native cultures."
  )
  print(explanation)

find_foreign_symbol()