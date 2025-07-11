import math

def calculate_kappa(tau_u, tau_v, phi, mu, rho):
    """
    Calculates the parameter kappa based on the model parameters.

    Args:
        tau_u (float): Time constant for the postsynaptic accumulator.
        tau_v (float): Time constant for the presynaptic accumulator.
        phi (float): Scaling constant for the presynaptic accumulator.
        mu (float): Mean firing rate of the synapses.
        rho (float): Offset constant in the Hebbian learning rule.

    Returns:
        float: The value of kappa.
    """
    kappa = - (tau_u + tau_v) * (phi * mu + rho) / phi
    return kappa

if __name__ == '__main__':
    # Example parameters for the model
    # Note: Rates should be in events per unit of time, consistent with time constants.
    # If time constants are in ms, rate should be in events/ms.
    tau_u = 20.0  # ms
    tau_v = 20.0  # ms
    phi = 1.0     # dimensionless scaling factor
    mu_hz = 5.0   # Hz
    mu = mu_hz / 1000.0  # Convert Hz to events/ms
    rho = -0.01   # Potentiation/depression offset (in same units as v_k)

    # Calculate kappa
    kappa_value = calculate_kappa(tau_u, tau_v, phi, mu, rho)

    # Print the definition and the calculation
    print("The definition of kappa is derived from the steady-state analysis of the synaptic plasticity model.")
    print("The formula for kappa is:")
    print("kappa = - (tau_u + tau_v) * (phi * mu + rho) / phi\n")
    
    print("Plugging in the example values:")
    print(f"  tau_u = {tau_u}")
    print(f"  tau_v = {tau_v}")
    print(f"  phi   = {phi}")
    print(f"  mu    = {mu} (from {mu_hz} Hz)")
    print(f"  rho   = {rho}\n")

    print("The calculation is:")
    print(f"kappa = - ({tau_u} + {tau_v}) * ({phi} * {mu} + {rho}) / {phi}")
    
    sum_tau = tau_u + tau_v
    pre_term = phi * mu
    full_numerator = -(sum_tau) * (pre_term + rho)
    
    print(f"kappa = - ({sum_tau}) * ({pre_term:.3f} + {rho}) / {phi}")
    print(f"kappa = - ({sum_tau}) * ({pre_term + rho:.3f}) / {phi}")
    print(f"kappa = ({full_numerator:.3f}) / {phi}")
    print(f"\nFinal calculated kappa: {kappa_value:.4f}")