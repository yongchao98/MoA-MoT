import math
from scipy.stats import ncx2

def solve_f_statistic_problem():
    """
    Calculates the minimum F-statistic for a given confidence level and max relative bias.
    """
    # Define problem parameters
    max_relative_bias = 0.10
    confidence_level = 0.95
    num_instruments = 1
    significance_level = 1 - confidence_level

    print("Here is the step-by-step derivation to find the minimum F-statistic:")
    print("-" * 70)

    # Step 1: Relate Relative Bias to the First-Stage F-statistic
    print("Step 1: Relate Relative Bias (b) to the expected First-Stage F-statistic (E[F]).")
    print("Using the weak-instrument asymptotic approximation from Staiger & Stock (1997),")
    print("the relative asymptotic bias is approximated by the equation: b ≈ 1 / E[F].\n")

    # Step 2: Set the condition
    print("Step 2: Apply the desired constraint on the relative bias.")
    print(f"We want the relative bias to be less than {max_relative_bias} (or {max_relative_bias*100}%).")
    print(f"The condition is: 1 / E[F] < {max_relative_bias}")
    expected_F = 1 / max_relative_bias
    print(f"Solving for E[F], we get: E[F] > {expected_F:.1f}.\n")

    # Step 3: Relate E[F] to the non-centrality parameter (μ²)
    print("Step 3: Relate E[F] to the model's concentration parameter (μ²).")
    print("The expected value of the F-statistic is related to its non-centrality parameter, μ²,")
    print("and the number of instruments, L, by the equation: E[F] ≈ 1 + μ²/L.")
    print(f"With L = {num_instruments} instrument, this becomes: E[F] ≈ 1 + μ².")
    mu_squared_null = expected_F - 1
    print(f"Combining with Step 2, our condition becomes 1 + μ² > {expected_F:.1f}, which implies μ² > {mu_squared_null:.1f}.\n")

    # Step 4: Formulate the hypothesis test
    print("Step 4: Formulate the hypothesis test.")
    print(f"We test the null hypothesis that the instrument is 'weak' (H0: μ² <= {mu_squared_null:.1f})")
    print(f"against the alternative that it is 'strong' (H1: μ² > {mu_squared_null:.1f}).")
    print(f"To have {confidence_level*100:.0f}% confidence that the bias is less than 10%, we need a statistical")
    print(f"test with a significance level (alpha) of {significance_level:.2f}.\n")

    # Step 5: Calculate the critical F-value
    print("Step 5: Calculate the critical F-value.")
    print("The critical value is the threshold for our observed F-statistic. We find it from the")
    print("distribution of the F-statistic under the boundary of the null hypothesis, where μ² = {0:.1f}.".format(mu_squared_null))
    print(f"For a large sample and L={num_instruments} instrument, the F-statistic is distributed as a non-central")
    print(f"chi-squared variable with df=L={num_instruments} and non-centrality parameter nc=μ²={mu_squared_null:.1f}.")
    print(f"We must find the critical value F_crit such that P(F > F_crit) = {significance_level:.2f}, which is the")
    print(f"{confidence_level*100:.0f}th percentile of this distribution.\n")

    # Perform the calculation
    df = num_instruments
    nc = mu_squared_null
    q = confidence_level
    critical_value = ncx2.ppf(q, df, nc)

    print("--- Final Calculation ---")
    print(f"The critical F-statistic is calculated as the {int(q*100)}th percentile of a non-central")
    print(f"chi-squared distribution with {df} degree of freedom and non-centrality parameter {nc:.1f}.")
    print(f"F_crit = ncx2.ppf(q={q}, df={df}, nc={nc:.1f})")
    print(f"The calculated critical F-statistic is: {critical_value:.4f}")

    final_answer = round(critical_value)
    print(f"\nRounding this to the nearest whole number, we get the final answer.")
    print(f"Final Answer: {final_answer}")
    return final_answer

if __name__ == '__main__':
    solve_f_statistic_problem()