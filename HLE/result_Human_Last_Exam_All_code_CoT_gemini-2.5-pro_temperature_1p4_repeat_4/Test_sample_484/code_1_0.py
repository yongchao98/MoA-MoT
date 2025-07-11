def display_final_equation():
    """
    This function prints the derived steady-state equation for synaptic efficacy.
    
    The equation describes the rate of change of synaptic weight w_i as a function
    of the presynaptic activity v_i and a shared postsynaptic activity u.
    
    Variables in the equation:
    - tau_w: Time constant for synaptic efficacy.
    - w_i: Efficacy of synapse i.
    - beta: Parameter for LTP strength.
    - u: Shared postsynaptic accumulator (u = sum(w_j * v_j)).
    - rho: A constant combining model parameters (rho = (alpha - beta)*(1 - eta)/beta).
    - phi: Efficiency of MMP9-induced cleavage.
    - v_i: Presynaptic accumulator for synapse i (i.e., its firing rate).
    """
    
    # Define the symbols for the equation as strings
    tau_w = "τ_w"
    w_i_dot = "dw_i/dt"
    beta = "β"
    rho = "ρ"
    phi = "φ"
    u = "u"
    v_i = "v_i"
    
    # Construct the final equation string
    # We are printing each symbol/term as requested.
    # The final equation is: τ_w * dw_i/dt = β * u * (1 + ρ / (1 + φ * v_i))
    
    term1 = f"{tau_w} * {w_i_dot}"
    term2 = f"{beta} * {u}"
    term3 = "1"
    term4 = f"{rho}"
    term5 = f"1 + {phi} * {v_i}"

    final_equation = f"{term1} = {term2} * ({term3} + {term4} / ({term5}))"
    
    print("The derived expression for the change in synaptic efficacy is:")
    print(final_equation)
    
    print("\nWhere the variables are defined as:")
    print(f"v_i = x_i (presynaptic rate)")
    print(f"u = Σ_j(w_j * v_j) (postsynaptic accumulator)")
    print(f"And the constant ρ is defined as:")
    print("ρ = (α - β) * (1 - η) / β")

display_final_equation()