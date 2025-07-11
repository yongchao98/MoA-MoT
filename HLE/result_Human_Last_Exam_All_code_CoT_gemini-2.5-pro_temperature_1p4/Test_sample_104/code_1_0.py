import scipy.stats

def find_minimum_f_statistic():
    """
    Calculates the minimum F-statistic for a TSLS model with one instrument
    to have 95% confidence that the relative asymptotic bias is less than 10%.
    """
    
    # --- Parameters ---
    relative_bias_threshold = 0.10
    confidence_level = 0.95
    num_instruments = 1
    
    print("This script calculates the minimum required F-statistic based on the following parameters:")
    print(f"- Number of instruments (K): {num_instruments}")
    print(f"- Desired maximal relative bias: {relative_bias_threshold:.2f} (10%)")
    print(f"- Desired confidence level: {confidence_level:.2f} (95%)")
    print("-" * 50)

    # --- Step 1: Define the hypothesis in terms of the expected F-statistic (E[F]) ---
    # The relative asymptotic bias is approximated by 1 / E[F].
    # We test the null hypothesis H0: relative_bias >= 0.10 vs. H1: relative_bias < 0.10.
    # This is equivalent to H0: E[F] <= 1/0.10 vs. H1: E[F] > 1/0.10.
    # We use the boundary of the null hypothesis for the calculation.
    expected_F_H0 = 1 / relative_bias_threshold
    print("Step 1: Determine the expected F-statistic at the null hypothesis boundary.")
    print(f"The equation is: E[F] = 1 / relative_bias = 1 / {relative_bias_threshold:.2f} = {expected_F_H0:.0f}")
    print("-" * 50)
    
    # --- Step 2: Calculate the non-centrality parameter (lambda) ---
    # The F-statistic ~ (1/K) * Chi2(K, lambda), and E[F] = 1 + lambda / K.
    # From this, we derive lambda = (E[F] - 1) * K.
    non_centrality_param = (expected_F_H0 - 1) * num_instruments
    print("Step 2: Calculate the non-centrality parameter (lambda) of the chi-squared distribution.")
    print(f"The equation is: lambda = (E[F] - 1) * K = ({expected_F_H0:.0f} - 1) * {num_instruments} = {non_centrality_param:.0f}")
    print("-" * 50)

    # --- Step 3: Find the critical value from the distribution ---
    # The critical value is the 95th percentile of the non-central chi-squared
    # distribution with df = K and ncp = lambda.
    # For K=1, the distribution is simply Chi2(df=1, ncp=9).
    df = num_instruments
    critical_value = scipy.stats.ncx2.ppf(confidence_level, df, non_centrality_param)
    print("Step 3: Find the critical F-value, which is the 95th percentile of this distribution.")
    print(f"Calculated critical value (F-statistic): {critical_value:.4f}")
    print("-" * 50)
    
    # --- Step 4: Round to the nearest whole number ---
    final_answer = round(critical_value)
    print("Final Step: Round the result to the nearest whole number.")
    print(f"The minimum F-statistic required is {final_answer}.")
    
    return final_answer

if __name__ == '__main__':
    find_minimum_f_statistic()