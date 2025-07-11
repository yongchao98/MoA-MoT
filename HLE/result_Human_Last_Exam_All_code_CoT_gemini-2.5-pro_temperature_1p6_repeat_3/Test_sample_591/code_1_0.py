import sys

def get_kappa_definition():
  """
  Prints the definition of the parameter kappa (κ) based on the provided model of synaptic plasticity.
  """
  # Define the variables symbolically for the explanation.
  # These are the parameters from the model provided by the user.
  rho_symbol = "ρ"  # Offset constant in the Hebbian equation
  w_symbol = "w"    # Assumed stable synaptic efficacy
  nu_symbol = "ν"   # Density of synapses on the dendrite
  phi_symbol = "ϕ"  # Scaling constant for the presynaptic accumulator

  # The definition of kappa (κ) relates these fundamental model parameters.
  # It essentially captures the relative strength of the offset-driven (depressing)
  # part of the learning rule compared to the direct presynaptic drive.

  # Print the final equation for kappa
  print("The definition of κ in the expression for critical correlation is:")
  print(f"κ = -({rho_symbol} * {w_symbol} * {nu_symbol}) / {phi_symbol}")

  # Explain each variable in the equation
  print("\nWhere:")
  print(f"  {rho_symbol} (rho): is the offset constant in the Hebbian learning rule for synaptic efficacy.")
  print(f"  {w_symbol} (w): is the assumed equal and stable synaptic efficacy.")
  print(f"  {nu_symbol} (nu): is the linear density of synapses along the dendrite.")
  print(f"  {phi_symbol} (phi): is the scaling constant for the presynaptic accumulator.")

# Execute the function to display the answer.
get_kappa_definition()
