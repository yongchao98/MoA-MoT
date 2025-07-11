def print_kappa_definition():
    """
    This function prints the definition of the variable kappa (κ)
    based on the steady-state analysis of the provided synaptic plasticity model.
    """

    # The derived mathematical expression for kappa
    kappa_formula = "κ = - ( (ρ + φ*μ*τ_v) * (τ_u + τ_v) ) / (φ*τ_v)"

    print("The definition of κ in the expression for critical correlation is:")
    print(kappa_formula)

    print("\nWhere each variable in the equation represents:")
    print("ρ (rho): The offset constant in the Hebbian learning rule.")
    print("φ (phi): The scaling constant for the presynaptic accumulator.")
    print("μ (mu): The mean firing rate of the synapses.")
    print("τ_v (tau_v): The time constant of the presynaptic accumulator.")
    print("τ_u (tau_u): The time constant of the postsynaptic accumulator.")

if __name__ == '__main__':
    print_kappa_definition()