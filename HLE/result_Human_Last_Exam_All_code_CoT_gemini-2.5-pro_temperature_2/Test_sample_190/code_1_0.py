import sys
# Set a higher recursion limit for display purposes if needed, though not required for this calculation.
sys.setrecursionlimit(2000)

def solve_transience_problem():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    The solution is presented step-by-step.
    """

    print("To determine the condition for transience, we analyze the drift and the second moment of the increments.")
    print("-" * 70)

    # Step 1: Calculate the expected drift mu_k
    print("Step 1: Calculate the expected drift mu_k = E[X_{n+1} - k | X_n=k].")
    print("The transitions from state k are to k-2, k-1, k+1, k+2.")
    print("The increments are -2, -1, 1, 2.")
    print("Probabilities are P_k,k-2=1/4, P_k,k+2=1/4, P_k,k-1=1/4-c/k, P_k,k+1=1/4+c/k.")
    
    print("\nmu_k = (-2) * P_{k,k-2} + (-1) * P_{k,k-1} + (1) * P_{k,k+1} + (2) * P_{k,k+2}")
    print("mu_k = (-2) * (1/4) + (-1) * (1/4 - c/k) + (1) * (1/4 + c/k) + (2) * (1/4)")
    print("mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2")
    print("mu_k = 2*c / k")
    print("-" * 70)

    # Step 2: Extract the coefficient mu from the asymptotic form mu_k = mu/k
    print("Step 2: From the asymptotic form of the drift (mu/k), we identify the coefficient μ.")
    mu_from_criterion = "2c"
    print(f"The drift mu_k has the form μ/k, where μ = {mu_from_criterion}.")
    print("-" * 70)

    # Step 3: Calculate the second moment of the increment, sigma_sq_k
    print("Step 3: Calculate the second moment of the increment, sigma_sq_k = E[(X_{n+1}-k)^2 | X_n=k].")
    print("\nsigma_sq_k = (-2)^2*P_{k,k-2} + (-1)^2*P_{k,k-1} + (1)^2*P_{k,k+1} + (2)^2*P_{k,k+2}")
    print("sigma_sq_k = 4 * (1/4) + 1 * (1/4 - c/k) + 1 * (1/4 + c/k) + 4 * (1/4)")
    print("sigma_sq_k = 1 + (1/4 - c/k + 1/4 + c/k) + 1")
    print("sigma_sq_k = 1 + 1/2 + 1")
    second_moment_val_num = 5
    second_moment_val_den = 2
    print(f"sigma_sq_k = {second_moment_val_num}/{second_moment_val_den}")
    print("-" * 70)

    # Step 4: Find the limit sigma^2
    print("Step 4: The limit of the second moment as k -> infinity gives us σ².")
    sigma_sq = f"{second_moment_val_num}/{second_moment_val_den}"
    print(f"σ² = lim_{{k->inf}} sigma_sq_k = {sigma_sq}")
    print("-" * 70)
    
    # Step 5: Apply the transience criterion
    print("Step 5: Apply the criterion for transience, which is 2*μ > σ².")
    lhs_coeff = 2
    
    print(f"Substituting μ = {mu_from_criterion} and σ² = {sigma_sq}, we get the inequality:")
    # Represents 2 * (2c) > 5/2
    final_lhs_coeff = 4
    print(f"{lhs_coeff} * ({mu_from_criterion}) > {sigma_sq}")
    print(f"This simplifies to: {final_lhs_coeff}c > {sigma_sq}")
    print("-" * 70)

    # Step 6: Solve for c and find the infimum
    print("Step 6: Solve the inequality for c to find the set where the chain is transient.")
    # c > (5/2) / 4 = 5/8
    result_num = second_moment_val_num
    result_den = second_moment_val_den * final_lhs_coeff
    
    print(f"c > ({sigma_sq}) / {final_lhs_coeff}")
    print(f"c > {result_num}/{result_den}")
    print("\nThe set of c values for which the chain is transient is ({}/{result_den}, infinity).".format(result_num, result_den=result_den))
    
    infimum = f"{result_num}/{result_den}"
    print(f"The infimum of this set is {infimum}.")
    
# Execute the function to display the solution.
solve_transience_problem()

print("<<<5/8>>>")