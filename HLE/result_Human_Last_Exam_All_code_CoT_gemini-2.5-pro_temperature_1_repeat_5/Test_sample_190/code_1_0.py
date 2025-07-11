import sys

def solve_markov_transience():
    """
    Calculates the infimum of c for a Markov chain to be transient.
    The steps are printed to explain the derivation.
    """

    print("Step 1: Calculate the asymptotic mean drift, mu_k = E[X_{n+1} - k | X_n = k]")
    print("The transition probabilities for large k are:")
    print("P(k, k+2) = 1/4")
    print("P(k, k-2) = 1/4")
    print("P(k, k+1) = 1/4 + c/k")
    print("P(k, k-1) = 1/4 - c/k")
    print("\nThe possible jumps are +2, -2, +1, -1.")
    print("mu_k = (+2)*P(k,k+2) + (-2)*P(k,k-2) + (+1)*P(k,k+1) + (-1)*P(k,k-1)")
    print("mu_k = 2*(1/4) - 2*(1/4) + 1*(1/4 + c/k) - 1*(1/4 - c/k)")
    print("mu_k = (1/4 + c/k) - (1/4 - c/k) = 2c/k")
    print("The drift is of the form A/k, where A = 2c.\n")
    A_coeff = 2  # A = 2*c

    print("Step 2: Calculate the asymptotic variance of the jumps, sigma^2")
    print("First, we calculate the second moment of the jumps, M2_k = E[(X_{n+1} - k)^2 | X_n = k].")
    print("M2_k = (+2)^2*P(k,k+2) + (-2)^2*P(k,k-2) + (+1)^2*P(k,k+1) + (-1)^2*P(k,k-1)")
    print("M2_k = 4*(1/4) + 4*(1/4) + 1*(1/4 + c/k) + 1*(1/4 - c/k)")
    second_moment = 4 * (1/4) + 4 * (1/4) + (1/4) + (1/4)
    print(f"M2_k = 1 + 1 + 1/4 + 1/4 = {second_moment}")
    print("The variance sigma_k^2 = M2_k - (mu_k)^2 = 2.5 - (2c/k)^2.")
    print("As k -> infinity, the term (2c/k)^2 -> 0.")
    sigma_sq = second_moment
    print(f"So, the asymptotic variance sigma^2 = {sigma_sq}.\n")

    print("Step 3: Apply the transience criterion")
    print("The criterion for transience is 2*A > sigma^2.")
    lhs_coeff = 2 * A_coeff
    print(f"Substituting A = {A_coeff}c and sigma^2 = {sigma_sq}, we get the inequality:")
    print(f"{2} * ({A_coeff}*c) > {sigma_sq}")
    print(f"{lhs_coeff}*c > {sigma_sq}\n")
    
    print("Step 4: Solve for c")
    # 4c > 2.5
    # c > 2.5 / 4
    infimum_c = sigma_sq / lhs_coeff
    print(f"c > {sigma_sq} / {lhs_coeff}")
    print(f"c > {infimum_c}\n")

    print("Step 5: Find the infimum")
    print("The set of values for c for which the chain is transient is (0.625, infinity).")
    print("The infimum of this set is the lower bound.")
    print(f"The infimum value is {infimum_c}.")

# Execute the function
solve_markov_transience()
# The final answer in the required format.
# Using sys.stdout.write to avoid extra newlines from print().
sys.stdout.write("<<<0.625>>>")