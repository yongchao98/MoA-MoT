import sympy

def get_kappa_definition():
    """
    This function provides the definition of the parameter kappa from the plasticity model.
    """
    # Define the symbols used in the model equations provided by the user.
    # tau_v: time constant of the presynaptic accumulator v_k
    # tau_u: time constant of the postsynaptic accumulator u_k
    tau_v = sympy.Symbol('τ_v')
    tau_u = sympy.Symbol('τ_u')
    kappa = sympy.Symbol('κ')

    # In the original paper, the expression for critical correlation c* is derived as:
    # c* = ((tau_v / tau_u) * S - 1) / (S - 1)
    # The user's prompt gives the formula:
    # c* = (kappa * S - 1) / (S - 1)
    # By comparing these two, we can define kappa.
    
    definition = sympy.Eq(kappa, tau_v / tau_u)
    
    print("The parameter 'kappa' (κ) is defined as the ratio of the presynaptic and postsynaptic time constants.")
    print("\nFrom the model equations:")
    print(f"τ_v is the time constant of the presynaptic accumulator v_k.")
    print(f"τ_u is the time constant of the postsynaptic accumulator u_k.")
    print("\nThe definition is:")
    print(definition)

if __name__ == '__main__':
    get_kappa_definition()