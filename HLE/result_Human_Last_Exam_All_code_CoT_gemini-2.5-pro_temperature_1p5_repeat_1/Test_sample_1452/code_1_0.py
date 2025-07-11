def get_born_equation():
  """
  This function prints the name and famous equation for Max Born.
  """
  name = "Max Born"
  
  # The Born rule states that the probability density (rho) of a particle
  # is proportional to the square of the magnitude of its wave function (Psi).
  # The equation contains the number 2.
  equation_latex = "\\rho = |\\Psi|^2"
  
  print(f"Name: {name}")
  print(f"Equation: \\( {equation_latex} \\)")

get_born_equation()