import math

def solve_imitation_learning_bound():
    """
    This function derives and prints the tightest upper bound for the performance
    difference J(pi^*) - J(pi_hat) based on the problem statement.
    """

    # Symbolic representation of the problem's parameters
    H = "H"
    A_size = "|A|"
    lam = "lambda"
    pi_star = "pi*"
    pi_hat = "pi_hat"

    print("Derivation of the upper bound for J(pi^*) - J(pi_hat):")
    print("-" * 50)

    # Step 1: Start with the standard performance bound for Behavioral Cloning.
    print(f"Step 1: The standard performance difference bound in imitation learning for a finite horizon '{H}' is:")
    print(f"J({pi_star}) - J({pi_hat}) <= {H}^2 * E[TV({pi_star}, {pi_hat})]")
    print("where E[TV(...)] is the average Total Variation distance between the policies under the expert's state distribution.")
    print("-" * 50)

    # Step 2: Interpret the given information on the population total variation risk.
    print("Step 2: We are given the risk inequality:")
    print(f"T({pi_hat}, {pi_star}) <= {A_size} * (1 - exp(-{lam}))")
    print("Given that the average TV distance is at most 1, we interpret the 'population total variation risk' T as a scaled version of the average TV distance:")
    print(f"T({pi_hat}, {pi_star}) = {A_size} * E[TV({pi_star}, {pi_hat})]")
    print("\nThis implies:")
    print(f"{A_size} * E[TV({pi_star}, {pi_hat})] <= {A_size} * (1 - exp(-{lam}))")
    print("=> E[TV({pi_star}, {pi_hat})] <= 1 - exp(-{lam})")
    print("-" * 50)
    
    # Step 3: Combine the results to find the final upper bound.
    print("Step 3: Substituting the bound on the average TV distance into the performance difference formula gives the final result.")

    # Final Equation components, as requested by the prompt to output "each number in the final equation"
    term1 = f"{H}^2"
    term2 = f"(1 - exp(-{lam}))"
    
    print("\n--- Final Upper Bound ---")
    print(f"The tightest upper bound on J({pi_star}) - J({pi_hat}) is the product of two terms:")
    print(f"1. Horizon-dependent term: {term1}")
    print(f"2. Algorithm-dependent risk term: {term2}")
    
    final_bound = f"{term1} * {term2}"
    print(f"\nFinal Equation: J({pi_star}) - J({pi_hat}) <= {final_bound}")

solve_imitation_learning_bound()