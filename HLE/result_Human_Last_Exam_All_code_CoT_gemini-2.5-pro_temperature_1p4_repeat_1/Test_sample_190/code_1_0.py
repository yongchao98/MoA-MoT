def find_infimum_c():
    """
    This script calculates the infimum of c for which the given Markov chain is transient.
    It follows a step-by-step analytical approach based on drift analysis.
    """
    
    print("This problem requires finding the threshold for 'c' that makes a Markov chain transient.")
    print("We will use a drift-based criterion (a form of Lamperti's Test).\n")

    # --- Step 1: Calculate the Mean Drift (mu_k) ---
    print("Step 1: Calculate the mean drift, mu_k = E[X_{n+1} - k | X_n = k].")
    print("The jumps from state k are {-2, -1, 1, 2}.")
    print("The corresponding probabilities are {1/4, 1/4 - c/k, 1/4 + c/k, 1/4}.")
    print("mu_k = (-2) * (1/4) + (2) * (1/4) + (-1) * (1/4 - c/k) + (1) * (1/4 + c/k)")
    print("mu_k = -1/2 + 1/2 - 1/4 + c/k + 1/4 + c/k")
    print("After simplification, the mean drift is: mu_k = 2c/k\n")

    # --- Step 2: Calculate the Second Moment of Jumps (sigma_k^2) ---
    print("Step 2: Calculate the second moment of jumps, sigma_k^2 = E[(X_{n+1} - k)^2 | X_n = k].")
    print("sigma_k^2 = (-2)^2 * (1/4) + (2)^2 * (1/4) + (-1)^2 * (1/4 - c/k) + (1)^2 * (1/4 + c/k)")
    print("sigma_k^2 = (4 * 1/4) + (4 * 1/4) + (1/4 - c/k) + (1/4 + c/k)")
    print("sigma_k^2 = 1 + 1 + 1/2")
    sigma_k_sq = 2.5
    print(f"After simplification, the second moment of the jumps is: sigma_k^2 = {sigma_k_sq}\n")
    
    # --- Step 3: Apply the Transience Criterion ---
    print("Step 3: Apply the transience criterion.")
    print("A random walk is transient if for large k, the condition 2*k*mu_k / sigma_k^2 > 1 holds.")
    print("Let's evaluate the expression: 2 * k * (2c/k) / (5/2)")
    print("This simplifies to: 4c / (5/2), which is 8c/5.")
    print("So, the condition for transience becomes the inequality: 8c/5 > 1\n")
    
    # --- Step 4: Solve the Inequality for c ---
    print("Step 4: Solve the inequality for c.")
    print("The final inequality to solve is:")
    print(8, "* c /", 5, ">", 1)
    print("Multiplying by 5 gives: 8 * c > 5")
    print("Dividing by 8 gives: c > 5/8\n")

    # --- Step 5: Determine the Infimum ---
    print("Step 5: Determine the infimum.")
    print("The chain is transient for all c in the set (5/8, infinity).")
    infimum = 5 / 8
    print(f"The infimum (greatest lower bound) of this set is 5/8, which is {infimum}.")

if __name__ == '__main__':
    find_infimum_c()
    print("<<<0.625>>>")
