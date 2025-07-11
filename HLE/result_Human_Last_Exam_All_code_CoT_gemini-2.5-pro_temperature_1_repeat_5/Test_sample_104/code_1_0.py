import math
from scipy.stats import ncx2

def solve_f_statistic_problem():
    """
    Calculates the minimum F-statistic for a given confidence level and relative bias.
    """
    # --- Problem Parameters ---
    # We have one instrumental variable (K2) and one endogenous variable.
    num_instruments = 1

    # We want 95% confidence that the relative asymptotic bias is less than 10%.
    max_relative_bias = 0.10
    confidence_level = 0.95
    significance_level = 1 - confidence_level

    # --- Step 1: Find the required minimum concentration parameter (lambda) ---
    # The relative bias is approximated by: Relative Bias ≈ 1 / E[F]
    # The expected F-statistic is approximated by: E[F] ≈ 1 + lambda / num_instruments
    # We want: 1 / (1 + lambda / num_instruments) < max_relative_bias
    # Solving for lambda:
    # 1 / max_relative_bias < 1 + lambda / num_instruments
    # lambda > num_instruments * (1 / max_relative_bias - 1)
    min_expected_f = 1 / max_relative_bias
    min_lambda = num_instruments * (min_expected_f - 1)

    print("--- Step 1: Determine the required instrument strength (Concentration Parameter λ) ---")
    print(f"The desired maximum relative bias is {max_relative_bias:.0%}.")
    print(f"This requires the expected F-statistic, E[F], to be greater than 1 / {max_relative_bias:.2f} = {min_expected_f:.1f}.")
    print(f"Since E[F] ≈ 1 + λ/K₂, with K₂ = {num_instruments}, we have the equation: 1 + λ/{num_instruments} > {min_expected_f:.1f}")
    print(f"Solving for λ gives the condition: λ > {min_lambda:.1f}.")
    print("-" * 60)


    # --- Step 2: Calculate the critical F-statistic ---
    # We test H₀: λ ≤ 9 against H₁: λ > 9.
    # The critical value is the 95th percentile of the F-statistic's distribution
    # under the worst-case null, where λ = 9.
    # Asymptotically, this F-statistic follows a non-central chi-squared distribution
    # with K₂ degrees of freedom and non-centrality parameter λ.
    df = num_instruments
    nc = min_lambda
    critical_f_statistic = ncx2.ppf(confidence_level, df=df, nc=nc)
    final_answer = round(critical_f_statistic)

    print("--- Step 2: Calculate the Critical F-Statistic for the Hypothesis Test ---")
    print(f"We need to be {confidence_level:.0%} confident that λ > {min_lambda:.1f}.")
    print(f"This requires our observed F-statistic to be greater than the critical value from a test where H₀: λ = {min_lambda:.1f}.")
    print(f"The critical value is the {confidence_level:.0%} percentile of the non-central chi-squared distribution with:")
    print(f"  - Degrees of Freedom (K₂) = {df}")
    print(f"  - Non-centrality Parameter (λ) = {nc:.1f}")
    print("\nThe calculated critical F-statistic is: {:.2f}".format(critical_f_statistic))
    print("-" * 60)

    # --- Final Answer ---
    print("--- Final Answer ---")
    print("The minimum F-statistic on the instrument such that you have 95% confidence that the")
    print("relative asymptotic bias is less than 10% is the calculated value rounded to the nearest whole number.")
    print(f"\nNearest whole number = {final_answer}")


if __name__ == "__main__":
    solve_f_statistic_problem()