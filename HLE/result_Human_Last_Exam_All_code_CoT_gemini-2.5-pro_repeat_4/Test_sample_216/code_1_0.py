import math

def solve_imitation_learning_bound():
    """
    This function derives and prints the tightest upper bound on the performance difference
    J(pi^*) - J(pi_hat) based on the provided information.
    """

    # --- Symbolic Representation ---
    # We use string representations for the mathematical symbols.
    H = "H"
    abs_A = "|A|"
    lambda_param = "lambda"
    pi_star = "pi^*"
    pi_hat = "pi_hat"
    J_pi_star = f"J({pi_star})"
    J_pi_hat = f"J({pi_hat})"
    R_max = "R_max"

    # --- Step-by-step Derivation ---
    print("Deriving the tightest upper bound for J(pi^*) - J(pi_hat):")
    print("-" * 60)

    # Step 1: State the performance difference bound
    print("Step 1: The Performance Difference Bound")
    print("A standard result in imitation learning theory bounds the difference in performance by the average one-step policy error.")
    print("This error measures how differently the learned policy behaves compared to the expert, averaged over the states the expert visits.")
    print("The tightest worst-case bound for an algorithm like Behavioral Cloning is quadratic in the episode horizon H:")
    print(f"  {J_pi_star} - {J_pi_hat} <= {H}^2 * {R_max} * E_s[d_TV({pi_star}(.|s), {pi_hat}(.|s))]")
    print(f"Here, {R_max} is the maximum single-step reward, and E_s[d_TV(...)] is the population total variation risk.\n")

    # Step 2: Use the given information
    print("Step 2: Incorporate the Given Risk Bound")
    print("The problem states that for this particular algorithm, the population TV risk is bounded:")
    risk_bound_expr = f"{abs_A} * (1 - exp(-{lambda_param}))"
    print(f"  E_s[d_TV({pi_star}(.|s), {pi_hat}(.|s))] <= {risk_bound_expr}\n")

    # Step 3: Combine and Conclude
    print("Step 3: Combine Bounds and Finalize")
    print("By substituting the given risk bound into the performance difference inequality, we get the final upper bound.")
    print(f"It is common practice to assume rewards are normalized to [0, 1], so we set {R_max} = 1.")
    
    # --- Final Equation ---
    final_bound_expr = f"{H}^2 * {abs_A} * (1 - exp(-{lambda_param}))"
    final_equation = f"{J_pi_star} - {J_pi_hat} <= {final_bound_expr}"

    print("\nThe tightest upper bound is:")
    print(final_equation)
    
    print("\n--- Components of the Final Equation ---")
    print(f"The performance gap '{J_pi_star} - {J_pi_hat}' is less than or equal to the product of:")
    print(f"1. Horizon Factor: {H}^2 (The number '2' is the exponent)")
    print(f"2. Action Space Factor: {abs_A}")
    print(f"3. Algorithm Performance Factor: (1 - exp(-{lambda_param})) (The number '1' is in this term)")
    print("-" * 60)

# Execute the function to print the solution
solve_imitation_learning_bound()