def print_derived_equation():
    """
    Prints the derived expression for the change in synaptic efficacy.
    """
    # The derived expression for the dynamics of synaptic efficacy w_i is:
    # τ_w * (dw_i/dt) = u_i * β * (v_i - ρ) / (1 + v_i)
    #
    # where:
    # w_i: synaptic efficacy
    # τ_w: time constant for synaptic efficacy
    # u_i: the postsynaptic accumulator (steady-state postsynaptic calcium Y)
    # v_i: the presynaptic accumulator (steady-state MMP9 level M_i)
    # β: relative strength of LTP
    # ρ: a constant threshold defined as ρ = - (α(1-η) + βη) / β
    #    α: relative strength of LTD
    #    η: baseline BDNF to proBDNF ratio parameter

    equation_str = "τ_w * dw_i/dt = u_i * β * (v_i - ρ) / (1 + v_i)"
    
    print("The derived expression for the dynamics of synaptic efficacy is:")
    print(equation_str)
    
    print("\nWhere the variables and constants are defined as:")
    print("τ_w: Time constant for synaptic efficacy")
    print("w_i: Synaptic efficacy")
    print("u_i: Postsynaptic accumulator (proportional to sum of all weighted inputs)")
    print("v_i: Presynaptic accumulator (proportional to presynaptic firing rate)")
    print("β: Relative strength of LTP")
    print("ρ: Constant plasticity threshold, defined in terms of the original model parameters α, β, and η")
    print("Each term in the final equation is:")
    print(f"Term 1 (LHS): τ_w * dw_i/dt")
    print(f"Term 2 (RHS): u_i")
    print(f"Term 3 (RHS): β")
    print(f"Term 4 (RHS, numerator): v_i - ρ")
    print(f"Term 5 (RHS, denominator): 1 + v_i")


print_derived_equation()
<<<τ_w*dw_i/dt = u_i*β*(v_i - ρ)/(1 + v_i)>>>