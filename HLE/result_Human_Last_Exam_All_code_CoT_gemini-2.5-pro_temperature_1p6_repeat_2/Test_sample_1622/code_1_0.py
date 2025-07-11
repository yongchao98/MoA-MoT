import math

def get_p_n_formula():
  """
  This function prints the formula for P(n).
  L represents ln(n).
  The formula is derived to make the relative error of the approximation
  for Q(n) be of the order O(((ln n)/n)^4).
  """
  
  # The formula is P(n) = term1 + term2
  # term1 is the coefficient for 1/n^2
  # term2 is the coefficient for 1/n^3
  
  numerator1 = "3*L**2 - 2*L + 2"
  denominator1 = "24*n**2"
  
  numerator2 = "L**3 - 2*L**2 + 2*L"
  denominator2 = "48*n**3"

  print("The formula for P(n) is:")
  print(f"P(n) = ({numerator1}) / ({denominator1}) + ({numerator2}) / ({denominator2})")

  # To fulfill the requirement of printing each number in the equation.
  print("\nHere is the equation with its numerical components explicitly shown:")
  print("P(n) = (", 3, "*L**2 - ", 2, "*L + ", 2, ") / (", 24, "*n**2) + (L**3 - ", 2, "*L**2 + ", 2, "*L) / (", 48, "*n**3)", sep='')

# Execute the function to display the formula.
get_p_n_formula()
