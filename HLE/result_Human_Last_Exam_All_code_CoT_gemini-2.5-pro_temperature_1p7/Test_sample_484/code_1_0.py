def print_final_equation():
    """
    This function prints the derived equation for the steady-state weight dynamics.
    The derived equation is: τ_w * dw_i/dt = u_i * (ρ + β*v_i) / (1 + v_i)

    Variables are defined as:
    - w_i: synaptic efficacy
    - v_i = ɸ * ν_i (steady-state presynaptic accumulator M_i)
    - u_i = Σ_j w_j * ν_j (steady-state postsynaptic accumulator Y)
    - ρ = α(1-η) + βη (a constant combining plasticity parameters)
    """
    
    # Print each symbol of the final equation to stdout.
    print("τ_w", "*", "dw_i/dt", "=", "u_i", "*", "(ρ", "+", "β", "*", "v_i)", "/", "(1", "+", "v_i)")

print_final_equation()