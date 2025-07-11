def get_kappa_definition():
  """
  This function constructs and prints the mathematical definition of kappa
  based on the derivation from the dendritic plasticity model.
  """

  # Define the parameters as symbolic strings
  kappa = "kappa"
  tau_u = "tau_u"
  tau_v = "tau_v"
  rho = "rho"
  phi = "phi"

  # Construct the equation string. The equation clearly shows each parameter.
  equation = f"{kappa} = -({tau_u} + {tau_v}) * (1 + {rho} / {phi})"

  # Print the final equation
  print(equation)

get_kappa_definition()