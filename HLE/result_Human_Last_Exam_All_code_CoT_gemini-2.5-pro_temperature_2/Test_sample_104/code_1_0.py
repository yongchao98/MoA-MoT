import numpy as np
from scipy.stats import ncf

def find_min_f_statistic():
    """
    Calculates the minimum F-statistic for a TSLS model with one instrument
    to have 95% confidence that the relative asymptotic bias is less than 10%.
    """
    # 1. Define problem parameters
    max_relative_bias = 0.10
    confidence_level = 0.95
    num_instruments = 1
    
    significance_level = 1 - confidence_level

    # 2. Determine the required expectation of the F-statistic
    # Relative Bias is approx. 1 / E[F]. We want Relative Bias < 0.10.
    # This implies 1 / E[F] < 0.10, which means E[F] > 10.
    # The boundary case for our null hypothesis (instrument is weak) is E[F] = 10.
    required_expected_f = 1 / max_relative_bias
    
    # 3. Determine the parameters of the non-central F-distribution
    # We are testing H0: E[F] <= 10. We use the boundary E[F] = 10.
    # The F-statistic follows a non-central F(dfn, dfd, nc) distribution.
    dfn = num_instruments
    # Assume a large sample for the denominator degrees of freedom
    dfd = 10000 
    
    # The non-centrality parameter (nc, or lambda) is related to E[F] by:
    # E[F] approx= 1 + nc / L for large dfd. (Here L = dfn).
    # nc = L * (E[F] - 1)
    non_centrality_param = dfn * (required_expected_f - 1)

    # 4. Calculate the critical value of the F-statistic
    # This is the (1 - alpha) quantile of the non-central F-distribution.
    # We reject H0 if our observed F-statistic is greater than this value.
    min_f_stat = ncf.ppf(confidence_level, dfn=dfn, dfd=dfd, nc=non_centrality_param)
    
    # 5. Round the result to the nearest whole number as requested
    rounded_f_stat = round(min_f_stat)

    # 6. Print the explanation and the result
    print("This script calculates the minimum first-stage F-statistic needed to ensure, with 95% confidence, that the TSLS relative bias is below 10%.")
    print("\n--- Calculation Steps ---")
    print(f"1. Maximum acceptable relative bias (d): {max_relative_bias}")
    print(f"2. Required minimum expectation of F-statistic: E[F] > 1/d = {required_expected_f:.2f}")
    print(f"3. Null hypothesis to reject: H0: E[F] <= {required_expected_f:.2f}")
    print(f"4. Confidence level = {confidence_level*100}%, Significance level (alpha) = {significance_level:.2f}")
    print("\n--- Distribution Parameters ---")
    print(f"The F-statistic follows a non-central F-distribution, F(dfn, dfd, nc)")
    print(f"- Numerator degrees of freedom (dfn, num instruments): {dfn}")
    print(f"- Denominator degrees of freedom (dfd, large sample): {dfd}")
    print(f"- Non-centrality parameter (nc): {non_centrality_param:.2f}, derived from E[F] and dfn")
    print("\n--- Final Calculation ---")
    print(f"To satisfy the condition, the observed F-statistic must be greater than the critical value,")
    print(f"which is the {confidence_level*100:.0f}th percentile of the F({dfn}, {dfd}, nc={non_centrality_param:.2f}) distribution.")
    print(f"\nCalculated Minimum F-statistic = {min_f_stat:.2f}")
    print(f"Rounding to the nearest whole number, the answer is: {rounded_f_stat}")


if __name__ == '__main__':
    find_min_f_statistic()