import scipy.stats
import numpy as np

def find_minimum_f_statistic():
    """
    Calculates the minimum F-statistic for a TSLS model with one instrument
    to have 95% confidence that the relative asymptotic bias is less than 10%.
    """

    # --- Parameters ---
    # The maximum tolerated relative bias for the TSLS estimator.
    target_relative_bias = 0.10
    # The confidence level required to ensure the bias is below the target.
    confidence_level = 0.95
    # The number of instruments in the model.
    num_instruments = 1

    # --- Derivation ---
    # The relative bias 'r' is approximated by the equation: r ≈ 1 / (1 + λ)
    # where λ is the non-centrality parameter of the F-statistic's distribution.
    # We solve for λ at the boundary where the relative bias is exactly the target.
    # Equation: λ = (1 / r) - 1
    lambda_boundary = (1 / target_relative_bias) - 1

    # The problem is now to find a critical F-statistic 'F_crit' such that we can
    # reject the null hypothesis H0: λ ≤ lambda_boundary with 95% confidence.
    # This F_crit is the 95th percentile of the F-statistic's distribution,
    # assuming λ is at its boundary value.
    # Asymptotically, for one instrument, the F-statistic follows a non-central
    # chi-squared distribution with df=1 and non-centrality parameter λ.
    df = num_instruments
    quantile = confidence_level

    # Calculate the critical value using SciPy's percentile point function (ppf).
    critical_f_statistic = scipy.stats.ncx2.ppf(quantile, df, lambda_boundary)
    final_answer = np.round(critical_f_statistic)

    # --- Output ---
    print("--- Calculation of Minimum First-Stage F-Statistic ---")
    print("\nParameters Used:")
    print(f"Target Relative Bias (r): {target_relative_bias}")
    print(f"Confidence Level: {confidence_level}")
    print(f"Number of Instruments (L): {num_instruments}")

    print("\nFinal Equation for Critical F-Statistic (F_crit):")
    print("F_crit = PPF(q, df, nc)")
    print("where PPF is the Percentile Point Function of the non-central chi-squared distribution.")
    print("\nNumbers in the Final Equation:")
    print(f"q (quantile): {quantile}")
    print(f"df (degrees of freedom): {df}")
    print(f"nc (non-centrality parameter, λ): {lambda_boundary}")

    print("\n--- Result ---")
    print(f"The calculated minimum F-statistic is: {critical_f_statistic:.4f}")
    print(f"The value rounded to the nearest whole number is: {final_answer:.0f}")

if __name__ == '__main__':
    find_minimum_f_statistic()
<<<17>>>