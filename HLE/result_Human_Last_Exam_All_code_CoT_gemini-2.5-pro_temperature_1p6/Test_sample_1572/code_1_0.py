def print_formula():
  """
  This function prints the final formula for P(n).
  The numbers in the formula are 1, 720, 1, and 5040.
  """
  numerator1 = 1
  denominator1 = 720
  
  numerator2 = 1
  denominator2 = 5040
  
  formula = f"P(n) = exp({numerator1}/({denominator1}*n**2) - {numerator2}/({denominator2}*n**4))"
  print(formula)

print_formula()