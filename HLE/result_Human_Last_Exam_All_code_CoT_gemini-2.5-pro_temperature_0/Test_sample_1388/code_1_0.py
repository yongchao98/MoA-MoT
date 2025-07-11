import math

def get_H_t_expression():
  """
  This function prints the explicit form of H(t) based on the mathematical derivation.
  The final form is H(t) = exp(h(t)/2).
  """
  
  # The final equation is H(t) = exp( (1/2) * h(t) )
  # The numbers in the equation are 1 and 2.
  numerator = 1
  denominator = 2
  
  print("The explicit form of H(t) is derived from the L2 energy estimate.")
  print("The final expression for H(t) is:")
  print(f"H(t) = exp(h(t) / {denominator})")
  # To explicitly show both numbers as requested
  print(f"   or H(t) = exp( ({numerator}/{denominator}) * h(t) )")

get_H_t_expression()