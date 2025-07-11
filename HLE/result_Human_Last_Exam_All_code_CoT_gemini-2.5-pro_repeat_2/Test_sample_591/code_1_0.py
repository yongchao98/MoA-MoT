import sys

def get_kappa_definition():
    """
    This function prints the definition of the constant kappa (κ) and explains its components.
    The definition is derived from the steady-state analysis of the provided dendritic plasticity model.
    """
    
    # The definition of kappa (κ)
    kappa_definition = "κ = -ρ / (φ * μ * τ_v)"

    # Explanation of the parameters in the formula, based on the provided text.
    # The prompt assumes equal firing rates μ_k = μ for all synapses.
    explanations = {
        "ρ (rho)": "The offset constant in the Hebbian learning rule for synaptic efficacy (w_k).",
        "φ (phi)": "The scaling constant for the presynaptic input in the presynaptic accumulator (v_k).",
        "μ (mu)": "The mean firing rate of the presynaptic neurons.",
        "τ_v (tau_v)": "The time constant for the presynaptic accumulator (v_k)."
    }

    print("The definition of κ (kappa) in the expression for critical correlation is:")
    print("-----------------------------------------------------------------------")
    print(kappa_definition)
    print("\nWhere the parameters are defined as:")
    for param, desc in explanations.items():
        print(f"- {param}: {desc}")
    print("-----------------------------------------------------------------------")

# Execute the function to print the definition
get_kappa_definition()
