def display_solution():
  """
  This function prints the derived formula for the ratio BM/MI.
  The side lengths a, b, and c are opposite to vertices A, B, and C, respectively.
  """
  
  # The formula is derived using the Incenter-Excenter Lemma and Ptolemy's Theorem.
  # The ratio is expressed in terms of the side lengths a, b, and c.
  
  formula_numerator = "a + c"
  formula_denominator = "b"

  print("The ratio BM/MI is expressed in terms of the side lengths a, b, and c as follows:")
  print(f"BM / MI = ({formula_numerator}) / {formula_denominator}")
  print("\nBreaking down the final equation:")
  print(f"Numerator: The sum of the lengths of side 'a' (BC) and side 'c' (AB).")
  print(f"Denominator: The length of side 'b' (AC).")

display_solution()