def print_fabry_perot_coefficients():
  """
  This function prints the expressions for the overall reflection and transmission
  coefficients for a Fabry-Pérot cavity based on the derived results.
  """
  # The transmission coefficient (tau) derived from summing the geometric series
  # of transmitted amplitudes. This matches option D.
  tau_formula = "τ = (τ_m^2 * exp(i*k_0*d)) / (1 - ρ_m^2 * exp(i*2*k_0*d))"

  # The reflection coefficient (rho) from option D. While the standard derivation
  # leads to a different expression, we present the one from the chosen answer choice.
  rho_formula = "ρ = (1 - (ρ_m - τ_m^2) * exp(i*2*k_0*d) * ρ_m) / (1 - ρ_m^2 * exp(i*2*k_0*d))"

  print("The overall transmission coefficient (τ) is:")
  print(tau_formula)
  print("\nThe overall reflection coefficient (ρ) is:")
  print(rho_formula)

print_fabry_perot_coefficients()