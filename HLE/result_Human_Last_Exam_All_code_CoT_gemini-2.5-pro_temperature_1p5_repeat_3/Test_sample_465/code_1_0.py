def solve_chromatic_roots_question():
  """
  Analyzes statements about chromatic and orbital chromatic roots and prints the final answer.

  - A: Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
  - B: Chromatic roots may not be real. (True)
  - C: Real chromatic roots may take on negative values. (True)
  - D: Real chromatic roots may take on non-integer values. (True)
  - E: Chromatic polynomials may have roots between 0 and 1. (False)

  The true statements are A, B, C, D.
  The final answer is the sorted string of these letters.
  """
  true_statements = ['A', 'B', 'C', 'D']
  
  # Sort the letters alphabetically
  true_statements.sort()
  
  # Join them into a single string
  answer = "".join(true_statements)
  
  print(answer)

solve_chromatic_roots_question()