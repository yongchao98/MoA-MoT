import sys

def get_kappa_definition():
  """
  This function prints the definition of the parameter kappa (κ)
  based on the provided biophysical model.
  """

  # Define the components of the formula as strings
  numerator = "2 * (ρ + μ*φ)"
  # The term τ_uv represents the effective time constant for the correlation between u and v
  tau_uv = "(τ_u * τ_v) / (τ_u + τ_v)"
  denominator = "τ_uv * φ"

  # Assemble the final expression for kappa
  kappa_expression = f"κ = - ({numerator}) / ({denominator})"

  # Print the definition clearly
  print("The parameter κ in the expression for the critical correlation c* is defined as follows:")
  print("\n" + kappa_expression)
  print("\nwhere the parameters are:")
  print("  ρ: the offset constant in the Hebbian equation.")
  print("  μ: the mean firing rate of the synapses.")
  print("  φ: the scaling constant for the presynaptic accumulator.")
  print("  τ_u, τ_v: the time constants for the postsynaptic and presynaptic accumulators.")
  print(f"  τ_uv: the effective correlation time constant, defined as {tau_uv}.")
  
  # The number '2' is explicitly part of the equation as requested.
  print(f"\nNote: The number in the equation is {2}.")

# Execute the function to display the answer
get_kappa_definition()