def print_kappa_definition():
    """
    This function prints the symbolic definition of the variable kappa (κ)
    based on the parameters of the dendritic plasticity model.
    """

    # The definition of kappa is derived from the steady-state analysis of the model.
    # It is expressed in terms of the following model parameters:
    # mu (μ): The mean firing rate of the presynaptic neurons.
    # rho (ρ): The offset constant in the Hebbian learning rule.
    # phi (φ): The scaling constant for the presynaptic accumulator.
    # tau_u (τ_u): The time constant for the postsynaptic accumulator.
    # tau_v (τ_v): The time constant for the presynaptic accumulator.

    kappa_definition = "- (mu + rho / phi) * (1 / tau_u + 1 / tau_v)"

    print("The definition of kappa (κ) in the expression for critical correlation c* is:")
    print(f"κ = {kappa_definition}")
    print("\nWhere the symbols represent:")
    print("μ: mean firing rate")
    print("ρ: Hebbian rule offset constant")
    print("φ: presynaptic accumulator scaling constant")
    print("τ_u: postsynaptic accumulator time constant")
    print("τ_v: presynaptic accumulator time constant")

if __name__ == "__main__":
    print_kappa_definition()