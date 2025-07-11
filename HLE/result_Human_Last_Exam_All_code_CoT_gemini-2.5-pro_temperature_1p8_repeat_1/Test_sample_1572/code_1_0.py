import sys

def provide_formula():
  """
  This function prints the derived formula for P(n).
  The formula is based on the asymptotic expansion of the ratio Q(n)/T(n).
  """
  
  # Coefficients for P(n) = 1 + p2*n^(-2) + p4*n^(-4)
  # where p_i are the same as a_i from the expansion of Q(n)/T(n).
  
  p2_num = 1
  p2_den = 720
  
  p4_num = -1433
  p4_den = 7257600

  # To follow the instruction "Remember in the final code you still need to output
  # each number in the final equation!", we construct the string with each number explicitly.
  
  term1 = "1"
  term2 = f"+ ({p2_num}/{p2_den})*n^(-2)"
  term3 = f" + ({p4_num}/{p4_den})*n^(-4)"
  
  # A more readable version might be:
  # P(n) = 1 + (1/720)*n^{-2} - (1433/7257600)*n^{-4}

  # Print the final formula for the user.
  print(f"P(n) = {term1} {term2} {term3}")

provide_formula()
