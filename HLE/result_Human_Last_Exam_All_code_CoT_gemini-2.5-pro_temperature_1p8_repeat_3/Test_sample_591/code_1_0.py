def get_kappa_definition():
  """
  This function prints the symbolic definition of the variable kappa.
  """
  
  # The definition of kappa as derived from the steady-state analysis
  # of the synaptic plasticity model.
  kappa_definition = "- (tau_u + tau_v) * (mu + rho / phi)"
  
  print("Based on the derivation from the model's steady-state condition, the definition of kappa is:")
  print(f"kappa = {kappa_definition}")

if __name__ == "__main__":
  get_kappa_definition()