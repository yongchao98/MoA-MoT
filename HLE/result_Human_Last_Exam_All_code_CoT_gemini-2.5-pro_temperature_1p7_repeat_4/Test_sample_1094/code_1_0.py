def print_ac_loss_formula():
  """
  This function prints the normalized AC loss formula for an elliptical superconductor.

  The problem asks for the loss per cycle per unit length (Q) in the standard form:
  2 * pi * Q / (mu_0 * Ic^2)
  where:
  - i = Im / Ic is the normalized current (Im = current amplitude, Ic = critical current).
  - The formula is valid for i < 1.

  Based on the critical state model (Norris's solution), the final expression is derived and printed.
  """

  # Define the variable and the terms of the equation as strings for clarity
  normalized_loss_lhs = "2 * pi * Q / (mu_0 * Ic^2)"
  
  # The formula contains numbers (2, 1, 1, 1, 2) which are all included.
  equation_term_1 = "(1 - i) * ln(1 - i)"
  equation_term_2 = "(1 + i) * ln(1 + i)"
  equation_term_3 = "i**2" # or i^2
  
  # Construct the final equation string
  final_equation_rhs = f"2 * [ {equation_term_1} + {equation_term_2} - {equation_term_3} ]"
  
  print(f"The normalized AC loss is given by the following function of i = Im/Ic:")
  print(f"{normalized_loss_lhs} = {final_equation_rhs}")

# Execute the function to print the formula
print_ac_loss_formula()