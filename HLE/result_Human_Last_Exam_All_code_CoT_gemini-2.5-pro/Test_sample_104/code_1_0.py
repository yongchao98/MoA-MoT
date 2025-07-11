import scipy.stats
import math

def solve_f_statistic_puzzle():
    """
    Calculates the minimum F-statistic for a TSLS estimation to ensure
    95% confidence that the relative asymptotic bias is less than 10%.
    """

    # --- Parameters based on the problem statement ---
    # The desired maximum relative asymptotic bias.
    relative_bias_max = 0.10

    # The confidence level required.
    confidence_level = 0.95
    significance_level = 1 - confidence_level

    # The number of instrumental variables.
    num_instruments = 1

    # --- Step 1: Relate Relative Bias to the F-statistic's concentration parameter (lambda) ---
    # The relative bias (RB) is approximated by RB ≈ k / λ, where k is the number of instruments
    # and λ is the concentration parameter (non-centrality parameter), which is E[F].
    # For k=1, RB ≈ 1 / λ.
    # We want RB <= 0.10, which implies 1 / λ <= 0.10, or λ >= 10.
    non_centrality_param = num_instruments / relative_bias_max

    # --- Step 2: Set up the hypothesis test ---
    # H0: The instrument is "weak" (λ < 10, or Relative Bias > 10%).
    # H1: The instrument is "strong enough" (λ >= 10, or Relative Bias <= 10%).
    # We need to find the critical F-statistic that allows us to reject H0 with 95% confidence.
    # This critical value is found at the boundary of the null hypothesis, i.e., where λ = 10.

    # --- Step 3: Calculate the critical value ---
    # The F-statistic (with one instrument) asymptotically follows a non-central chi-squared
    # distribution with `df = num_instruments` and non-centrality parameter `λ`.
    # We need the critical value 'c' such that P(F > c | λ=10) = 0.05.
    # This is the 95th percentile of the distribution.
    df = num_instruments
    f_critical_value = scipy.stats.ncx2.ppf(q=confidence_level, df=df, nc=non_centrality_param)

    # Round the final answer to the nearest whole number as requested.
    final_answer = round(f_critical_value)

    # --- Step 4: Print the explanation and results ---
    print("Plan:")
    print("1. Use the approximation that Relative Bias (RB) is inversely related to the expected F-statistic (λ), RB ≈ 1/λ.")
    print(f"2. The requirement RB <= {relative_bias_max} implies λ >= 1/{relative_bias_max}, so λ >= {non_centrality_param}.")
    print(f"3. Frame this as a hypothesis test to reject the null hypothesis that λ < {non_centrality_param} with {confidence_level:.0%} confidence.")
    print("4. Calculate the critical F-value from the non-central chi-squared distribution at the null hypothesis boundary.")
    print("-" * 20)
    print("Calculation:")
    print(f"The critical F-statistic is the {confidence_level:.0%}-percentile of a non-central chi-squared distribution with:")
    print(f"  - Degrees of freedom (df): {df} (for {num_instruments} instrument)")
    print(f"  - Non-centrality parameter (λ): {non_centrality_param}")
    print(f"\nCalculated F-statistic: {f_critical_value:.4f}")
    print(f"The nearest whole number is: {final_answer}")
    print("-" * 20)

    # Final equation format as requested
    print("Final Equation:")
    print(f"F_statistic ≈ ppf_ncx2(q={confidence_level}, df={df}, nc={int(non_centrality_param)})")
    print(f"{final_answer} ≈ ppf_ncx2(q={0.95}, df={1}, nc={10})")


solve_f_statistic_puzzle()
<<<16>>>