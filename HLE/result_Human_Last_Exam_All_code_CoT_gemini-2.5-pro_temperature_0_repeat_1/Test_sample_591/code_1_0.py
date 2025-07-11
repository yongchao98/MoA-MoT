def print_kappa_definition():
    """
    This function prints the mathematical definition of the parameter kappa (κ)
    from the provided model of dendritic plasticity.
    """
    # Define the symbols for the variables in the equation for clarity.
    # These are based on the parameters described in the model.
    kappa_symbol = "κ"
    tau_u_symbol = "τ_u"
    tau_v_symbol = "τ_v"
    w_symbol = "w"
    mu_symbol = "μ"

    # The definition of kappa is derived from the stability analysis of the model's
    # steady-state solution. It relates the system's time constants to its
    # steady-state firing rate and synaptic efficacy.

    # Print the formula in a structured way.
    print("The definition of the parameter κ (kappa) is:")
    print("")
    print(f"      {tau_u_symbol} + {tau_v_symbol}")
    print(f"{kappa_symbol}  =  ―――――――――――――――――――――")
    print(f"    2 * {tau_u_symbol} * {tau_v_symbol} * {w_symbol} * {mu_symbol}")
    print("")

    # Explain each symbol in the equation as requested.
    print("Where each symbol in the final equation represents:")
    print(f"  {kappa_symbol}:\tA parameter from the stability analysis relating the variance and mean of the postsynaptic potential.")
    print(f"  {tau_u_symbol}:\tThe time constant of the postsynaptic accumulator u_k.")
    print(f"  {tau_v_symbol}:\tThe time constant of the presynaptic accumulator v_k.")
    print(f"  {w_symbol}:\tThe uniform steady-state synaptic efficacy.")
    print(f"  {mu_symbol}:\tThe uniform mean firing rate for the synapses.")

# Execute the function to print the definition.
print_kappa_definition()