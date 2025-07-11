import scipy.stats

def solve_f_statistic():
    """
    Calculates the minimum F-statistic for a TSLS estimator under specific
    bias and confidence constraints.
    """
    # Number of instrumental variables
    L = 1
    
    # Desired maximal relative asymptotic bias
    b = 0.10
    
    # Confidence level and corresponding significance level (alpha)
    confidence_level = 0.95
    alpha = 1 - confidence_level

    # The null hypothesis of the weak instrument test is that the relative bias
    # is greater than or equal to our threshold 'b'. This corresponds to a
    # null hypothesis on the non-centrality parameter 'lambda'.
    # H0: bias >= b  <=>  1 / (1 + lambda / L) >= b
    # This simplifies to lambda <= L * ((1/b) - 1)
    # The critical value is found under the worst-case null, at the boundary.
    non_centrality_parameter = L * ((1 / b) - 1)

    # The F-statistic is asymptotically distributed as a non-central chi-squared
    # distribution with L degrees of freedom (df) and the calculated
    # non-centrality parameter (nc).
    # We need to find the critical value F_crit such that P(F > F_crit) = alpha.
    # This is the (1 - alpha) percentile of the distribution.
    f_critical_value = scipy.stats.ncx2.ppf(1 - alpha, df=L, nc=non_centrality_parameter)

    # Round the final answer to the nearest whole number as requested.
    final_answer = round(f_critical_value)

    print("To find the minimum F-statistic, we solve for the critical value in a weak instrument test.")
    print("The critical value, F_crit, is found using the percentile point function (ppf) of a non-central chi-squared distribution.")
    print("\nThe equation is of the form: F_crit = ppf(q, df, nc)")
    print("\nThe numbers in the final equation are:")
    print(f"q (quantile) = 1 - alpha = {1 - alpha}")
    print(f"df (degrees of freedom) = Number of Instruments (L) = {L}")
    print(f"nc (non-centrality parameter) = lambda = {non_centrality_parameter}")
    
    print(f"\nThe calculated F-statistic is: {f_critical_value:.2f}")
    print(f"The minimum F-statistic rounded to the nearest whole number is: {final_answer}")

solve_f_statistic()