def print_fabry_perot_coefficients():
  """
  This function prints the overall reflection and transmission coefficients
  for the Fabry-Pérot cavity based on the solved problem.
  """

  # Define the expressions as strings
  # tau represents the overall transmission coefficient
  tau_expression = "τ = (τ_m^2 * e^(i*k_0*d)) / (1 - ρ_m^2 * e^(i*2*k_0*d))"

  # rho represents the overall reflection coefficient
  rho_expression = "ρ = (1 - (ρ_m - τ_m^2) * e^(i*2*k_0*d) * ρ_m) / (1 - ρ_m^2 * e^(i*2*k_0*d))"

  # Print the final expressions
  print("The overall transmission coefficient τ is:")
  print(tau_expression)
  print("\nThe overall reflection coefficient ρ is:")
  print(rho_expression)

# Execute the function to display the result
print_fabry_perot_coefficients()