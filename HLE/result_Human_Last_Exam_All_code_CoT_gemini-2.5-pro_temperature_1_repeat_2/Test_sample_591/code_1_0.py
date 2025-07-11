def define_kappa():
    """
    This function prints the definition of the kappa (κ) parameter from the synaptic stability model.
    """
    
    # Define the components of the kappa formula.
    # These are parameters from the model described in the problem.
    rho = "ρ"  # The offset constant in the Hebbian learning rule.
    W_thr = "W_thr"  # The synaptic efficacy threshold for removal.
    phi = "φ"  # The scaling constant for the presynaptic accumulator.

    # Explain the context
    print("In the context of the provided synaptic plasticity model, the term κ (kappa) is a dimensionless constant.")
    print("It represents a ratio of parameters that determines the balance between synaptic potentiation and depression, relative to the structural stability threshold.")
    print("\nIts definition is:")
    
    # Print the final equation with each component clearly labeled
    print(f"κ = -{rho} / ({W_thr} * {phi})")

    # Explain each symbol in the equation
    print("\nWhere:")
    print(f"  {rho} (rho): The offset constant in the Hebbian learning rule for synaptic efficacy w_k. It typically has a negative value to allow for depression.")
    print(f"  {W_thr}: The fixed threshold for synaptic efficacy. If a synapse's efficacy w_k falls below W_thr, it is removed.")
    print(f"  {phi} (phi): The scaling constant for the presynaptic accumulator v_k, linking it to the presynaptic activity x_k(t).")

# Execute the function to display the definition.
define_kappa()