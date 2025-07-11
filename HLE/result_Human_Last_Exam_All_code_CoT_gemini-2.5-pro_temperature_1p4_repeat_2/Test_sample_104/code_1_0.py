import scipy.stats

def solve_f_statistic_problem():
    """
    Calculates the minimum F-statistic for 95% confidence that TSLS relative
    asymptotic bias is less than 10%.
    """
    # Define parameters based on the problem and common econometric approximations.
    # Desired maximal relative bias = 10% (0.10)
    relative_bias_threshold = 0.1
    # We have one instrumental variable.
    num_instruments = 1
    # We require 95% confidence.
    confidence_level = 0.95

    # Step 1: From the relative bias condition, find the required expected F-statistic.
    # Relative Bias ≈ 1 / E[F] < 0.1  =>  E[F] > 10
    required_expected_f = 1 / relative_bias_threshold

    # Step 2: From the E[F] condition, find the required non-centrality parameter (NCP or lambda).
    # For L=1 instrument, E[F] ≈ 1 + lambda.
    # 1 + lambda > 10  =>  lambda > 9
    lambda_threshold = required_expected_f - num_instruments

    # Step 3: Find the critical F-value for a hypothesis test.
    # We test H₀: lambda <= 9 vs. H₁: lambda > 9 at 5% significance.
    # The critical value is the 95th percentile of the F-statistic's distribution
    # when lambda = 9. Asymptotically, this is a non-central chi-squared distribution.
    df = num_instruments
    critical_f_statistic = scipy.stats.ncx2.ppf(confidence_level, df=df, nc=lambda_threshold)

    # Round the final result to the nearest whole number.
    rounded_f_statistic = round(critical_f_statistic)

    # --- Output the reasoning and the result ---
    print("This problem asks for the minimum first-stage F-statistic that provides 95% confidence")
    print("that the relative asymptotic bias of the TSLS estimator is less than 10%.")
    print("\nHere is the step-by-step derivation:")

    print("\nStep 1: Relate Relative Bias to the Expected F-Statistic")
    print(f"Using the approximation 'Relative Bias ≈ 1 / E[F]', the condition 'Relative Bias < {relative_bias_threshold}'")
    print(f"implies that we need E[F] > 1 / {relative_bias_threshold}, or E[F] > {required_expected_f}.")

    print("\nStep 2: Relate Expected F-Statistic to the Non-Centrality Parameter (λ)")
    print(f"For {num_instruments} instrument, the asymptotic approximation is E[F] ≈ {num_instruments} + λ.")
    print(f"So, we require {num_instruments} + λ > {required_expected_f}, which simplifies to λ > {lambda_threshold}.")

    print("\nStep 3: Find the Critical F-Value using a Hypothesis Test")
    print("We want to be 95% confident that λ > 9. This corresponds to a hypothesis test:")
    print(f"  H₀ (Null): The instrument is weak (λ ≤ {lambda_threshold})")
    print(f"  H₁ (Alternative): The instrument is strong enough (λ > {lambda_threshold})")
    print("The critical F-value for a test with 5% significance (95% confidence) is the value")
    print(f"that cuts off the top 5% of the distribution when λ = {lambda_threshold}.")

    print("\nStep 4: The Final Calculation")
    print("Asymptotically, the F-statistic follows a non-central chi-squared distribution with")
    print(f"{df} degree of freedom and a non-centrality parameter (ncp) of {lambda_threshold}.")
    print("The final equation for the critical value is:")
    print(f"F_critical = 95th_percentile(chi-squared(df={df}, ncp={lambda_threshold}))")
    print(f"F_critical ≈ {critical_f_statistic:.3f}")

    print("\nAnswer:")
    print(f"The minimum F-statistic, rounded to the nearest whole number, is {int(rounded_f_statistic)}.")

if __name__ == '__main__':
    solve_f_statistic_problem()