import sys

def solve():
    """
    Calculates the infimum of c for which the Markov chain is transient.
    """

    # Let the jump from state k be denoted by the random variable D_k.
    # The possible values for D_k are -2, -1, +1, +2.
    # The probabilities for these jumps for a large state k are:
    # P(D_k = -2) = 1/4
    # P(D_k = -1) = 1/4 - c/k
    # P(D_k = +1) = 1/4 + c/k
    # P(D_k = +2) = 1/4

    print("Step 1: Calculate the asymptotic drift parameter 'theta'.")
    print("The drift mu_k is the expected value of the jump D_k.")
    print("mu_k = E[D_k] = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)")
    print("mu_k = -0.5 - 0.25 + c/k + 0.25 + c/k + 0.5 = 2*c/k")
    print("For a random walk with drift mu_k ~ theta/k, we have theta = 2*c.")
    theta_coeff_in_c = 2
    print(f"The parameter 'theta' is {theta_coeff_in_c} * c.\n")


    print("Step 2: Calculate the asymptotic variance of the jump 'sigma^2'.")
    print("The variance V_k is E[D_k^2] for large k (since E[D_k] -> 0).")
    # For large k, mu_k -> 0, so Var(D_k) = E[D_k^2] - (E[D_k])^2 approx E[D_k^2].
    print("sigma^2 = E[D_k^2] = (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4)")
    print("sigma^2 = 4*(0.25) + 1*(0.25 - c/k) + 1*(0.25 + c/k) + 4*(0.25)")
    print("sigma^2 = 1 + 0.25 + 0.25 + 1 = 2.5")
    sigma_sq = 2.5
    print(f"The asymptotic variance 'sigma^2' is {sigma_sq}.\n")

    print("Step 3: Apply Lamperti's criterion for transience.")
    print("The criterion for transience is: theta > sigma^2 / 2")
    # Using 'f-string' to show the numbers being substituted into the formula.
    print(f"Substituting the expressions for theta and sigma^2 gives: {theta_coeff_in_c} * c > {sigma_sq} / 2\n")

    print("Step 4: Solve the inequality for c to find the transience condition.")
    # The final equation for the boundary of c.
    print(f"The inequality simplifies to: c > {sigma_sq} / (2 * {theta_coeff_in_c})")
    
    # Calculate the values
    numerator = sigma_sq
    denominator = 2 * theta_coeff_in_c
    critical_c = numerator / denominator
    
    # Output the final calculation step
    print(f"c > {numerator} / {denominator}")
    print(f"c > {critical_c}\n")

    print("The set of c for which the Markov chain is transient is ({}, infinity).".format(critical_c))
    print("The infimum of this set is the lower bound of the interval.")
    
    # Print the final numerical answer.
    print(f"\nThe infimum is {critical_c}.")


solve()