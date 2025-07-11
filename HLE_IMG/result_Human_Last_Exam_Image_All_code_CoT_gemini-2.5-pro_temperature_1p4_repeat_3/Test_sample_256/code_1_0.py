def solve_puzzle():
  """
  Analyzes the image of rock art to find an out-of-place symbol.
  """
  is_true = True
  symbol_description = "The letters 'NU' (or similar shapes resembling letters from the Latin alphabet)."
  explanation = (
      "The statement is True. In the upper right section of the pictograph panel, "
      "a set of symbols that strongly resemble the letters 'NU' from the Latin alphabet can be distinguished. "
      "The styles of the other paintings are characteristic of ancient Southwest Native American cultures, "
      "which did not use the Latin alphabet. Therefore, these letters are a symbol not related to these cultures, "
      "likely added at a much later date."
  )
  
  print(f"{is_true}. The symbol is: {symbol_description}")
  print("\nExplanation:")
  print(explanation)

solve_puzzle()