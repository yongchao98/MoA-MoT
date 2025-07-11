import scipy.stats as stats
import numpy as np

def find_minimum_f_statistic():
    """
    Calculates the minimum F-statistic for a TSLS model with one instrument
    to have 95% confidence that the relative asymptotic bias is less than 10%.
    """
    
    # Step 1: Define the parameters for the calculation.
    # The problem translates to finding a critical value for a hypothesis test.
    # The null hypothesis is that the non-centrality parameter (lambda) is 9.
    # We test this at a 5% significance level.

    # Degrees of freedom for the non-central chi-squared distribution (number of instruments)
    df = 1
    
    # Non-centrality parameter (lambda) at the boundary of the null hypothesis.
    # This is derived from the 10% relative bias condition: 1/(1+lambda) < 0.1 => lambda > 9
    nc = 9.0
    
    # We want to be 95% confident, which corresponds to the 95th percentile of the distribution.
    quantile = 0.95
    
    # Step 2: Calculate the critical value of the F-statistic.
    # This value is the 95th percentile of the non-central chi-squared distribution
    # with df=1 and non-centrality parameter nc=9.
    # We use the Percent Point Function (ppf), which is the inverse of the Cumulative Distribution Function (CDF).
    f_statistic_critical_value = stats.ncx2.ppf(quantile, df=df, nc=nc)
    
    # Step 3: Round the result to the nearest whole number as requested.
    final_answer = int(np.round(f_statistic_critical_value))

    # Step 4: Print the explanation and the final result.
    print("The problem is to find the minimum F-statistic 'F' that satisfies the condition:")
    print("95% confidence that Relative Bias < 10%")
    print("\nThis condition is equivalent to testing the null hypothesis that the instrument is 'weak' (non-centrality parameter lambda <= 9) at a 5% significance level.")
    print("The equation for the F-statistic critical value 'c' is:")
    print(f"P(F > c | df={df}, non-centrality={nc}) = {1-quantile:.2f}")
    print("\nThis critical value 'c' is the 95th percentile of the corresponding non-central chi-squared distribution.")
    print(f"The calculated F-statistic critical value is: {f_statistic_critical_value:.3f}")
    print(f"\nRounding to the nearest whole number, the minimum required F-statistic is: {final_answer}")

if __name__ == '__main__':
    find_minimum_f_statistic()