def find_probability_formula():
  """
  This function prints the derived formula for the probability P_m.
  The derivation steps are outlined above.
  The final formula expresses P_m in terms of m.
  """
  
  # The numerator of the probability formula P_m
  numerator = "(m + 1) * (m + 2)"
  
  # The denominator of the probability formula P_m
  # Breaking it down to show all the numbers in the final equation
  term1 = "2"
  term2 = "(2 * m + 1)"
  term3 = "(4 * m + 1)"
  denominator = f"{term1} * {term2} * {term3}"

  # Print the final formula for P_m
  print("The probability P_m is given by the formula:")
  print(f"P_m = {numerator} / ({denominator})")

find_probability_formula()