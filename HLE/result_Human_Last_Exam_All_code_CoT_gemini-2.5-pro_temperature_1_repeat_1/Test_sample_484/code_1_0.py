def derive_synaptic_plasticity_rule():
    """
    This function prints the results of a steady-state analysis on the
    provided biophysical model of a dendritic segment.
    """

    # --- Define symbols as strings for clear output ---
    tau_w_dot_wi = "τ_w * dw_i/dt"
    beta = "β"
    u_i_str = "u_i"
    v_i_str = "v_i"
    rho_str = "ρ"
    alpha = "α"
    eta = "η"
    w_i = "w_i"
    x_i = "<x_i>"
    phi = "ϕ"

    # --- Construct the final equation string ---
    final_equation = f"{tau_w_dot_wi} = ({beta} * {u_i_str} * ({v_i_str} - {rho_str})) / (1 + {v_i_str})"

    # --- Construct the definition strings ---
    u_i_def = f"{u_i_str} = Y = Σ_j({w_i.replace('_i', '_j')} * {x_i.replace('_i', '_j')})"
    v_i_def = f"{v_i_str} = M_i = {phi} * {x_i}"
    rho_def = f"{rho_str} = -({eta} + ({alpha}/{beta}) * (1 - {eta}))"

    # --- Print the results ---
    print("The steady-state analysis reduces the complex system to a simpler one described by the following variables:")
    print(f"  - {w_i}: Synaptic efficacy.")
    print(f"  - {v_i_str}: The presynaptic accumulator, defined as the steady-state level of MMP9: {v_i_def}")
    print(f"  - {u_i_str}: The postsynaptic accumulator, defined as the shared postsynaptic calcium level: {u_i_def}")
    print("    (Note: u_i is a global variable, identical for all synapses.)")
    print("-" * 30)

    print("The derived expression for the dynamics of synaptic efficacy is:")
    print(f"\n  {final_equation}\n")

    print("Each term in the final equation represents:")
    print(f"  - {tau_w_dot_wi}: The rate of change of synaptic efficacy.")
    print(f"  - {beta}: The relative strength of LTP.")
    print(f"  - {u_i_str}: The postsynaptic activity level.")
    print(f"  - {v_i_str}: The presynaptic activity level.")
    print(f"  - {rho_str}: The modification threshold for presynaptic activity.")
    print("-" * 30)

    print("The constant ρ, which acts as a threshold, is defined by the underlying parameters as:")
    print(f"\n  {rho_def}\n")
    print(f"where {alpha} is the strength of LTD and {eta} is related to the baseline BDNF/proBDNF ratio.")


if __name__ == '__main__':
    derive_synaptic_plasticity_rule()
    print("\n<<<τ_w*dw_i/dt = (β*u_i*(v_i - ρ))/(1 + v_i)>>>")
