import scipy.stats

def solve_f_statistic_problem():
    """
    Calculates the minimum F-statistic for a TSLS estimation to have 95% confidence
    that the relative asymptotic bias is less than 10%.
    """
    # --- Introduction ---
    print("This program calculates the required minimum F-statistic based on a formal hypothesis test.")
    print("The problem: In a TSLS regression with one instrumental and one endogenous variable, find the minimum")
    print("first-stage F-statistic to be 95% confident that the relative asymptotic bias is less than 10%.")
    
    # --- Step 1: Relate Bias to the Concentration Parameter (μ) ---
    print("\nStep 1: Define the relationship between relative bias and instrument strength.")
    print("A common approximation for the relative asymptotic bias (b) under 'local-to-zero' asymptotics is:")
    print("b = 1 / (1 + μ/k)")
    print("where μ is the concentration parameter and k is the number of instruments.")
    
    k = 1
    max_bias = 0.10
    print(f"\nIn our case, we have k = {k} instrument and we want the bias b < {max_bias}.")
    print(f"So, the equation is: 1 / (1 + μ / {k}) < {max_bias}")
    print("Solving for μ gives: 1 + μ > 10, which means μ > 9.")
    target_mu = 9
    print(f"Thus, the instrument is considered 'strong enough' if its concentration parameter μ is greater than {target_mu}.")

    # --- Step 2: Formulate the Hypothesis Test ---
    print("\nStep 2: Frame the problem as a hypothesis test.")
    print("We want to be 95% confident that μ > 9. This means we need to reject the null hypothesis that the instrument is weak.")
    print(f"Null Hypothesis (H0): μ <= {target_mu} (The instrument is weak or on the boundary)")
    print(f"Alternative Hypothesis (H1): μ > {target_mu} (The instrument is strong enough)")
    
    significance_level = 0.05
    confidence_level = 1 - significance_level
    print(f"We will test this at a {significance_level*100:.0f}% significance level (corresponding to {confidence_level*100:.0f}% confidence).")

    # --- Step 3: Use the F-Statistic and its Distribution ---
    print("\nStep 3: Find the critical value using the F-statistic's distribution.")
    print("The first-stage F-statistic is our test statistic. Asymptotically, it follows a non-central chi-squared distribution divided by the degrees of freedom.")
    print("Distribution: F ~ χ²(df, nc) / df")
    
    df = k # The degrees of freedom (df) is the number of instruments.
    print(f"For k={k} instrument, the degrees of freedom df = {df}.")
    print("The F-statistic's distribution simplifies to: F ~ χ²(1, nc).")
    
    # The non-centrality parameter (nc) under the "worst-case" null is our target μ.
    nc = target_mu
    print(f"The worst case under the null hypothesis is μ = {target_mu}, so we set the non-centrality parameter nc = {nc}.")
    
    # The critical value is the 95th percentile of this distribution.
    quantile = confidence_level
    print(f"The critical value is the {int(quantile*100)}th percentile of the χ²(df={df}, nc={nc}) distribution.")

    # --- Step 4: Calculate the Final Answer ---
    print("\nStep 4: Final calculation.")
    critical_value = scipy.stats.ncx2.ppf(quantile, df=df, nc=nc)
    
    print(f"The calculation is: ppf({quantile}, df={df}, nc={nc})")
    print(f"The required F-statistic must be greater than {critical_value:.3f}.")
    
    final_answer = round(critical_value)
    print(f"\nThe nearest whole number to this F-statistic is {final_answer}.")
    
    return final_answer

if __name__ == '__main__':
    answer = solve_f_statistic_problem()
    # The final answer is wrapped according to the required format.
    # print(f"<<<{answer}>>>")