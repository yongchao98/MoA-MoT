import math

def display_final_equation():
  """
  Prints the derived formula for the force per unit area on the conducting plane.
  The formula is derived from electromagnetic principles as outlined in the plan.
  """

  # Define symbolic parts of the equation for clarity
  term_mu0 = "μ₀"
  term_K0_sq = "K₀²"
  term_sin_sq = "sin²(ay)"
  term_cosh = "cosh(ad)"
  term_sinh = "sinh(ad)"
  term_mu_ratio = "μ₀/μ"
  
  # Construct the numerator and denominator strings
  numerator = f"{term_mu0} * {term_K0_sq} * {term_sin_sq}"
  denominator_base = f"{term_cosh} + {term_mu_ratio} * {term_sinh}"
  denominator = f"2 * ({denominator_base})²"
  
  # The force is a vector in the negative x direction
  direction = "î_x"

  print("The force per unit y-z area on the x=d interface is given by the equation:")
  print(f"      - ( {numerator} )")
  print(f"f/A = ----------------------------------------- * {direction}")
  print(f"        ( {denominator} )")

display_final_equation()
