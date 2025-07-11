def solve_parentheses_statements():
  """
  This function determines the truth value of six statements about parenthesis strings
  and prints the result as a single string.
  
  The analysis is based on a counterexample string family Sk = "(()()...())" with k inner pairs.
  - Statement 1: False. For Sk, LHS is Theta(k), RHS is O(1).
  - Statement 2: False. For Sk, LHS grows with k, RHS is 0.
  - Statement 3: False. For Sk, LHS is Theta(k), RHS is O(1).
  - Statement 4: True. For Sk, both LHS and RHS are Theta(k). No counterexample found.
  - Statement 5: True. For Sk, both LHS and RHS are Theta(k). No counterexample found.
  - Statement 6: True. For Sk, both LHS and RHS are Theta(k). No counterexample found.
  
  The key distinction arises from the behavior of the functions at D=1.
  """
  
  # The results for each statement (F for False, T for True)
  results = "FFFTTT"
  
  print(results)

solve_parentheses_statements()