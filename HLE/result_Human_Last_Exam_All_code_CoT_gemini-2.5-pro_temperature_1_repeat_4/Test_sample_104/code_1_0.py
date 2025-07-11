import scipy.stats
import math

def solve_f_statistic_threshold():
    """
    Calculates the minimum F-statistic for a TSLS estimation with one instrument
    to have 95% confidence that the relative asymptotic bias is less than 10%.
    """

    # --- Parameters based on the problem statement ---
    # Number of instrumental variables
    num_instruments = 1
    # Maximum desired relative bias (10%)
    max_relative_bias = 0.10
    # Required confidence level
    confidence_level = 0.95
    # Significance level (alpha) for the hypothesis test
    alpha = 1 - confidence_level

    # --- Step 1: Derive the non-centrality parameter (NCP) under the null hypothesis ---
    # The relative bias (RB) is approximated by RB ≈ 1 / (1 + lambda/K).
    # We test the null hypothesis that instruments are "weak", meaning the relative
    # bias is greater than or equal to our threshold.
    # H0: RB >= max_relative_bias  =>  1 / (1 + lambda/K) >= 0.10
    # Solving for lambda gives the boundary condition for the null hypothesis.
    # 1 / max_relative_bias <= 1 + lambda / K
    # (1 / max_relative_bias) - 1 <= lambda / K
    lambda_over_K_null = (1 / max_relative_bias) - 1
    ncp_null = lambda_over_K_null * num_instruments

    # --- Step 2: Calculate the critical value for the F-statistic ---
    # Asymptotically, for K=1 instrument, the F-statistic follows a non-central
    # chi-squared distribution with K degrees of freedom and non-centrality parameter lambda.
    # We need to find the critical value 'c' such that P(F > c | lambda = ncp_null) = alpha.
    # This is the (1-alpha) or 95th percentile of the specified distribution.
    degrees_of_freedom = num_instruments

    # Calculate the critical value using SciPy's non-central chi-squared
    # percentile point function (ppf).
    critical_value_f_stat = scipy.stats.ncx2.ppf(
        confidence_level,
        df=degrees_of_freedom,
        nc=ncp_null
    )

    # --- Step 3: Round the result to the nearest whole number ---
    final_answer = round(critical_value_f_stat)

    # --- Output the results and reasoning ---
    print("--- TSLS First-Stage F-Statistic Calculation ---")
    print("\nParameters:")
    print(f"Number of instruments (K): {num_instruments}")
    print(f"Maximum desired relative bias: {max_relative_bias:.2f} (or {max_relative_bias*100}%)")
    print(f"Required confidence level: {confidence_level:.2f} (or {confidence_level*100}%)")

    print("\nDerivation:")
    print(f"The condition on relative bias [1 / (1 + lambda/K) < {max_relative_bias}] implies that lambda must be greater than a certain value.")
    print(f"This value is (1 / {max_relative_bias}) - 1 = {lambda_over_K_null}.")
    print(f"We test the null hypothesis that instruments are weak: H0: lambda <= {ncp_null}.")

    print("\nCalculation:")
    print("The minimum F-statistic is the critical value of a test of H0 at a 5% significance level.")
    print("This corresponds to the 95th percentile of a non-central chi-squared distribution with:")
    print(f"  - Degrees of freedom (df): {degrees_of_freedom}")
    print(f"  - Non-centrality parameter (nc): {ncp_null}")
    print(f"Calculated critical value (the required F-statistic): {critical_value_f_stat:.4f}")

    print("\nFinal Answer:")
    print(f"The minimum F-statistic rounded to the nearest whole number is: {final_answer}")

if __name__ == '__main__':
    solve_f_statistic_threshold()