def solve_steady_state_pi0():
    """
    This function explains the derivation of the steady-state probability pi_0
    for the given birth-death process and prints the final formula.
    """

    # Define rho as a string for display purposes
    rho = "ρ"
    lambda_sym = "λ"
    mu_sym = "μ"

    print("Derivation for the steady-state probability pi_0:")
    print("-" * 50)

    # Step 1: Define the relationship between pi_k and pi_0
    print("Step 1: Express pi_k in terms of pi_0.")
    print(f"The birth rates are lambda_i = {lambda_sym}/(i+1) and death rates are mu_i = {mu_sym}.")
    print(f"The general relationship is pi_k = pi_0 * Product_{{i=0 to k-1}}(lambda_i / mu_{{i+1}}).")
    print(f"Substituting the rates gives pi_k = pi_0 * ({rho}^k / k!), where {rho} = {lambda_sym}/{mu_sym}.")
    print("-" * 50)

    # Step 2: Use the normalization condition
    print("Step 2: Use the normalization condition Sum(pi_k) for k=0 to infinity = 1.")
    print(f"Sum(pi_0 * ({rho}^k / k!)) = 1")
    print(f"pi_0 * Sum(({rho}^k / k!)) = 1")
    print("-" * 50)

    # Step 3: Recognize the Taylor series for e^rho
    print("Step 3: Recognize the summation as the Taylor series for e^x.")
    print(f"The sum Sum(({rho}^k / k!)) from k=0 to infinity is equal to e^{rho}.")
    print(f"So, the equation becomes: pi_0 * e^{rho} = 1.")
    print("-" * 50)

    # Step 4: Final solution for pi_0
    print("Step 4: Solve for pi_0.")
    final_equation_lhs = "π_0"
    final_equation_rhs = f"e^(-{rho})"
    print(f"The final expression for the steady-state probability is:")
    print(f"{final_equation_lhs} = {final_equation_rhs}")

solve_steady_state_pi0()