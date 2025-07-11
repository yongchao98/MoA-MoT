def print_final_equation():
    """
    Prints the derived simplified equation for synaptic efficacy dynamics.
    """
    
    # Define the symbols used in the equation for clarity
    tau_w = "τ_w"
    dw_dt = "dw_i/dt"
    beta = "β"
    u_i = "u_i"
    rho = "ρ"
    phi = "φ"
    v_i = "v_i"
    alpha = "α"
    eta = "η"

    # The main differential equation
    main_equation = f"{tau_w} * {dw_dt} = {beta} * {u_i} * (1 - {rho} / (1 + {phi} * {v_i}))"

    # The definition of the postsynaptic accumulator u_i
    u_definition = f"where {u_i} = Σ_j(w_j * v_j)"
    
    # The definition of the constant rho
    rho_definition = f"and {rho} = (1 - {alpha}/{beta}) * (1 - {eta})"

    print("The final simplified expression for the synaptic efficacy dynamics is:")
    print(main_equation)
    print("")
    print("In this equation, the variables and constants are defined as:")
    print("- w_i: The synaptic efficacy of synapse i.")
    print("- v_i: The presynaptic accumulator, representing the mean firing rate of synapse i.")
    print(f"- {u_i}: The postsynaptic accumulator, a shared variable representing the total weighted input.")
    print(f"  {u_definition}")
    print(f"- {rho}: A dimensionless constant that combines the effects of LTD, LTP, and BDNF/proBDNF baseline ratio.")
    print(f"  {rho_definition}")

print_final_equation()