def print_derived_equation():
    """
    This function prints the final derived equation for the synaptic weight dynamics
    after a steady-state analysis, along with the definitions of the simplified variables.
    """
    
    # The derived relationship for the dynamics of synaptic efficacy (w_i).
    # tau_w * dw_i/dt = u_i * ( (beta * v_i + rho) / (1 + v_i) )
    # Note: we use unicode for Greek letters for better readability in terminals that support it.
    tau = "\u03C4"  # tau
    beta = "\u03B2" # beta
    rho = "\u03C1"  # rho
    
    final_equation = f"{tau}_w * dw_i/dt = u_i * (({beta} * v_i + {rho}) / (1 + v_i))"
    
    print("After performing a steady-state analysis on the biophysical model, we derive the following simplified expression for the dynamics of synaptic efficacy, w_i:")
    print("\n" + "="*80)
    print(final_equation)
    print("="*80 + "\n")
    print("Where the variables are defined as:")
    print("  - v_i: The presynaptic accumulator, representing the steady-state level of presynaptic MMP9.")
    print("         v_i = M_i")
    print("\n  - u_i: The postsynaptic accumulator, representing the shared postsynaptic calcium level.")
    print("         u_i = Y")
    print("\n  - rho (\u03C1): A constant defined by a combination of the model's intrinsic parameters.")
    print("         \u03C1 = \u03B1(1 - \u03B7) + \u03B2\u03B7")
    
print_derived_equation()