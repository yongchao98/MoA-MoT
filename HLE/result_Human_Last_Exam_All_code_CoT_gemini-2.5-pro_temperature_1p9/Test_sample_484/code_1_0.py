def solve_synaptic_dynamics():
    """
    This function presents the result of a steady-state analysis on the provided
    biophysical model of a dendritic segment.

    The original system is reduced to a single equation describing the slow
    dynamics of the synaptic weight, w_i, in terms of a presynaptic
    accumulator, v_i, and a postsynaptic accumulator, u_i.

    Variable definitions for the final equation:
    - w_i: The synaptic efficacy (weight) of synapse i.
    - u_i: The postsynaptic accumulator, representing the steady-state
           postsynaptic calcium level, Y.
    - v_i: The presynaptic accumulator, representing the steady-state
           MMP9 level at synapse i, M_i.
    - rho: A constant defined as -alpha/beta, representing the relative
           strength of LTD over LTP.
    - beta: The parameter for the strength of LTP.
    - eta: The parameter determining the baseline BDNF/proBDNF ratio.
    """

    # The derived equation for the dynamics of synaptic efficacy w_i
    # The term (rho*(1-eta) - eta) acts as a plasticity threshold for v_i.
    final_equation = "τ_w * dw_i/dt = β * u_i * (v_i - (ρ*(1-η) - η)) / (1 + v_i)"

    print("Derived Equation for Synaptic Efficacy Dynamics:\n")
    print(final_equation)
    print("\nWhere each term is defined as follows:")
    print("----------------------------------------")
    print("     w_i: Synaptic efficacy")
    print("   dw_i/dt: Rate of change of synaptic efficacy")
    print("     τ_w: Time constant for synaptic efficacy")
    print("     u_i: Postsynaptic accumulator (steady-state Calcium Y)")
    print("     v_i: Presynaptic accumulator (steady-state MMP9 M_i)")
    print("       β: LTP strength parameter")
    print("       ρ: LTD/LTP ratio constant (-α/β)")
    print("       η: Baseline BDNF to proBDNF ratio parameter")
    print("----------------------------------------")
    print("\nFull equation with symbols printed out explicitly:")
    tau_w = "τ_w"
    dw_dt = "dw_i/dt"
    beta = "β"
    u_i = "u_i"
    v_i = "v_i"
    rho = "ρ"
    eta = "η"
    one = "1"
    print(f"{tau_w} * {dw_dt} = {beta} * {u_i} * ({v_i} - ({rho} * ({one} - {eta}) - {eta})) / ({one} + {v_i})")


solve_synaptic_dynamics()

# <<<τ_w * dw_i/dt = β * u_i * (v_i - (ρ*(1-η) - η)) / (1 + v_i)>>>