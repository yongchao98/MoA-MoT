def display_final_equation():
    """
    Prints the derived equation for synaptic efficacy dynamics.
    """
    
    # Define the components of the equation as strings
    tau_w_dot_w_i = "τ_w * dw_i/dt"
    beta_u_i = "β * u_i"
    rho = "ρ"
    v_i = "v_i"
    
    # Construct the final equation string
    final_equation = f"{tau_w_dot_w_i} = {beta_u_i} * (1 - {rho} / (1 + {v_i}))"
    
    # Define the expression for rho
    rho_definition = f"where ρ = (1 - η) * (β - α) / β"
    
    # Print the results
    print("The derived expression for the dynamics of synaptic efficacy is:")
    print(final_equation)
    print("\n" + rho_definition)
    
    # Explain the variables
    print("\nWhere:")
    print("w_i: Synaptic efficacy")
    print("u_i: Postsynaptic accumulator (proportional to total weighted postsynaptic input)")
    print("v_i: Presynaptic accumulator (proportional to presynaptic firing rate)")
    print("α, β, η: Model parameters")


if __name__ == "__main__":
    display_final_equation()