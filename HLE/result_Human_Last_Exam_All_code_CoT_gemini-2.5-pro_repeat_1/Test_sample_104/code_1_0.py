from scipy.stats import ncx2

def find_minimum_f_statistic():
    """
    Calculates the minimum F-statistic for 95% confidence that the relative 
    asymptotic bias is less than 10% in a TSLS estimation with one instrument.
    """
    
    # --- Step 1: Define parameters ---
    relative_bias_threshold = 0.10
    confidence_level = 0.95
    alpha = 1 - confidence_level
    num_instruments = 1

    # --- Step 2: Calculate the Non-Centrality Parameter (NCP) for the bias threshold ---
    # The relationship is: Relative Bias ≈ 1 / (1 + NCP)
    # We solve for NCP at the bias threshold of 10% (0.10).
    ncp = (1 / relative_bias_threshold) - 1
    
    print("Goal: Find the F-statistic (F) for 95% confidence that Relative Bias < 10%.")
    print("\nStep 1: Relate Relative Bias to the Non-Centrality Parameter (NCP).")
    print(f"Relative Bias ≈ 1 / (1 + NCP)")
    print(f"{relative_bias_threshold} = 1 / (1 + NCP)")
    print(f"Solving for NCP gives: NCP = (1 / {relative_bias_threshold}) - 1 = {ncp}")

    # --- Step 3: Find the critical F-value ---
    # We need to find the F-statistic that allows us to reject the null hypothesis
    # that NCP <= 9 (i.e., bias >= 10%) with 95% confidence.
    # This is the 95th percentile of the distribution under the null hypothesis,
    # which is a non-central chi-squared distribution with df=1 and nc=9.
    degrees_of_freedom = num_instruments
    
    # The ppf (Percent Point Function) is the inverse of the CDF.
    # ppf(0.95) gives the value below which 95% of the distribution lies.
    critical_f_value = ncx2.ppf(confidence_level, df=degrees_of_freedom, nc=ncp)
    
    print("\nStep 2: Find the critical F-value from the corresponding distribution.")
    print("We need the F-value that puts 5% in the upper tail of the distribution")
    print("defined by the NCP boundary. This is the 95th percentile.")
    print(f"Distribution: Non-central chi-squared(df={degrees_of_freedom}, nc={ncp})")
    print(f"F_critical = ppf({confidence_level}, df={degrees_of_freedom}, nc={ncp})")
    print(f"F_critical ≈ {critical_f_value:.4f}")

    # --- Step 4: Round to the nearest whole number and report ---
    final_answer = round(critical_f_value)
    
    print(f"\nStep 3: Round the result to the nearest whole number.")
    print(f"The minimum required F-statistic is {final_answer}.")
    
    return final_answer

# Run the function to get the final answer.
# The function prints the logic and the final answer will be captured below.
final_result = find_minimum_f_statistic()
# The required output format is just the final answer wrapped in <<<>>>
# print(f'<<<{final_result}>>>')