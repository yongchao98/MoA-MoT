def print_final_equation():
    """
    Prints the derived steady-state equation for synaptic efficacy.
    """
    # Define the variables and constants for clarity in the printout
    tau_w = "tau_w"
    w_i_dot = "dw_i/dt"
    u_i = "u_i"
    v_i = "v_i"
    beta = "beta"
    rho = "rho"
    
    # The final derived equation
    final_equation = f"{tau_w} * {w_i_dot} = {u_i} * ({beta} + {rho} / (1 + {v_i}))"
    
    # Definition of the constant rho
    rho_definition = "where rho = (1 - eta) * (alpha - beta)"

    # Print the results
    print("The derived expression for the dynamics of synaptic efficacy is:")
    print(final_equation)
    print("\nWith the constant rho defined as:")
    print(rho_definition)
    print("\nAnd the variables defined as:")
    print(f"{u_i}: postsynaptic accumulator (Y)")
    print(f"{v_i}: presynaptic accumulator (M_i)")


# Execute the function to display the answer
print_final_equation()