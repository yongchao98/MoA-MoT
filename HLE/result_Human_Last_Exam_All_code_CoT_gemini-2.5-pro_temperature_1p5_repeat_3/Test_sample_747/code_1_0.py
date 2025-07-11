def solve():
  """
  This function provides the solution to the user's request.
  It identifies the incorrect lines in the provided proof and determines the correctness of the claim.
  The final answer is printed to the console.
  """
  # List of line numbers containing wrong statements.
  incorrect_lines = [3, 7]
  # The claim itself is wrong.
  claim_correctness = "Wrong"
  
  # The final output format is a list followed by a word.
  print(f"{incorrect_lines} {claim_correctness}")

solve()