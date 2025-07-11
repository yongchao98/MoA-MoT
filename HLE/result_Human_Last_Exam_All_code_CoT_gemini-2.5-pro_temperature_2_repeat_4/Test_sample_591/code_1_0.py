import math

def calculate_kappa(rho, tau_u, phi, tau_v):
  """
  Calculates the parameter kappa (κ).

  κ represents the ratio of the postsynaptic offset term to the
  presynaptic scaling term in the weight dynamics.
  """
  return - (rho * tau_u) / (phi * tau_v)

def calculate_c_star(kappa, S):
  """
  Calculates the critical correlation c* for synaptic stabilization.
  """
  if S == 1:
    return float('inf') # Avoid division by zero
  return (kappa * S - 1) / (S - 1)

def main():
  """
  Main function to define parameters and print the solution.
  """
  # --- Definition of Kappa ---
  # By deriving the critical correlation from the model equations,
  # kappa (κ) is found to be defined as:
  print("The parameter kappa (κ) is a dimensionless quantity defined as the ratio of key model parameters:")
  print("κ = - (ρ * τ_u) / (φ * τ_v)\n")

  # --- Example Calculation ---
  # Let's use some plausible example values for the parameters.
  rho = -0.5      # Synaptic offset constant
  tau_u = 40.0    # Postsynaptic accumulator time constant (ms)
  phi = 1.0       # Presynaptic accumulator scaling constant
  tau_v = 20.0    # Presynaptic accumulator time constant (ms)
  S = 10.0        # Effective number of neighboring synapses spatially coupled

  print("Let's calculate κ using the following example parameters:")
  print(f"ρ (rho)     = {rho}")
  print(f"τ_u (tau_u) = {tau_u}")
  print(f"φ (phi)     = {phi}")
  print(f"τ_v (tau_v) = {tau_v}\n")

  # Calculate kappa
  kappa = calculate_kappa(rho, tau_u, phi, tau_v)

  # Print the kappa calculation
  print("The calculation for κ is:")
  print(f"κ = - ({rho} * {tau_u}) / ({phi} * {tau_v})")
  print(f"κ = {kappa}\n")

  # Calculate c_star
  c_star = calculate_c_star(kappa, S)

  # Print the c_star calculation
  print("Now, we use this κ to find the critical correlation c* with S = {S}:".format(S=S))
  print(f"c* = (κ * S - 1) / (S - 1)")
  print(f"c* = ({kappa} * {S} - 1) / ({S} - 1)")
  print(f"c* = ({kappa * S} - 1) / ({S - 1})")
  print(f"c* = {kappa * S - 1} / {S - 1}")
  print(f"c* = {c_star}")

if __name__ == "__main__":
  main()