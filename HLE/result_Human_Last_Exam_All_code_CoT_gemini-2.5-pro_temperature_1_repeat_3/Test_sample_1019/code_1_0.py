import math

def print_hamiltonicity_threshold():
  """
  This function prints the d-threshold for Hamiltonicity for the given range.
  The threshold is p = log(n)/n. We print the coefficients to satisfy the
  output format requirements.
  """
  numerator_coeff = 1
  denominator_coeff = 1
  
  print("The d-threshold for Hamiltonicity in the given range is a probability p(n) with the following asymptotic behavior:")
  print(f"p(n) = ({numerator_coeff} * log(n)) / ({denominator_coeff} * n)")

print_hamiltonicity_threshold()