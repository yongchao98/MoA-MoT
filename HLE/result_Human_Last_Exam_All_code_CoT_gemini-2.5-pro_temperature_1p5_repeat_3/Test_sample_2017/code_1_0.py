def solve_greek_puzzle():
  """
  This function identifies non-classical words in two Greek passages,
  provides their classical equivalents, and the periods they belong to.
  """
  # Analysis for Passage 1
  w1 = "εἶμαι"
  c1 = "εἰμί"
  p1 = "KoineDemotic"

  # Analysis for Passage 2
  w2 = "ὄρνιθα"
  c2 = "ὄρνιν"
  p2 = "Homeric"

  # Combine the results into the specified format
  result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
  
  print(result)

solve_greek_puzzle()