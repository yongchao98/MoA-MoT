import math

def display_hamiltonicity_threshold():
  """
  This function constructs and prints the symbolic formula for the d-threshold for Hamiltonicity.
  """
  
  # Define the symbolic variables and numbers in the equation
  n_variable = "n"
  eta_variable = "eta"
  number_two = 2
  
  # Construct the expression for the minimum degree 'd'
  d_expression = f"({n_variable} / {number_two} - {eta_variable})"
  
  # Construct the numerator and denominator of the threshold probability formula
  numerator = f"ln({d_expression})"
  denominator = d_expression
  
  # Construct the full equation as a string
  final_equation = f"p({n_variable}) = {numerator} / {denominator}"
  
  print("The d-threshold for Hamiltonicity in this range is given by the formula:")
  print(final_equation)

if __name__ == "__main__":
  display_hamiltonicity_threshold()