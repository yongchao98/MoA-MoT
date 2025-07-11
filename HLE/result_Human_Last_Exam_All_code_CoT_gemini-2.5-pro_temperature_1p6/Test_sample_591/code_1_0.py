def get_kappa_definition():
  """
  This function prints the definition of the variable kappa.
  """
  # In the provided model, kappa is a parameter derived from other constants.
  # It is defined by the ratio of two key constants from the model equations.

  # 'rho' is the offset constant from the synaptic efficacy (w_k) update rule.
  # 'phi' is the scaling constant from the presynaptic accumulator (v_k) equation.
  
  kappa_definition = "-rho / phi"

  print("The definition of kappa is given by the following equation:")
  print(f"kappa = {kappa_definition}")

get_kappa_definition()