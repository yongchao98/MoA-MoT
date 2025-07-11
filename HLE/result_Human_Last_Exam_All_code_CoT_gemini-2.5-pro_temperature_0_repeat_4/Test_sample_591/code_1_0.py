def define_kappa():
    """
    Prints the definition of the parameter kappa (κ) from the model of synaptic plasticity.
    """
    definition = "κ = - (ρ / (φ * μ * τ_v)) * ((τ_u + τ_v) / τ_u)"
    
    print("The definition of κ (kappa) is:")
    print(definition)
    print("\nWhere the parameters are:")
    print("ρ (rho): The offset constant in the Hebbian learning rule.")
    print("φ (phi): The scaling constant for the presynaptic accumulator.")
    print("μ (mu): The mean firing rate of the presynaptic neurons.")
    print("τ_v (tau_v): The time constant of the presynaptic accumulator.")
    print("τ_u (tau_u): The time constant of the postsynaptic accumulator.")

if __name__ == "__main__":
    define_kappa()