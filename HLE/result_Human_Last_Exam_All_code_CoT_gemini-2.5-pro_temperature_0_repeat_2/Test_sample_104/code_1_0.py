import scipy.stats
import math

def calculate_minimum_f_statistic():
    """
    Calculates the minimum F-statistic for a given relative bias and confidence level.
    """
    # --- Step 1: Define parameters ---
    relative_bias_threshold = 0.10
    confidence_level = 0.95
    num_instruments = 1
    significance_level = 1 - confidence_level

    # --- Step 2: Use the approximation to find the required population F-statistic ---
    # Relative Bias ≈ 1 / E[F]
    # We want Relative Bias < 0.10, so E[F] > 1 / 0.10
    required_population_F = 1 / relative_bias_threshold
    
    print("Step 1: Determine the required population F-statistic (E[F])")
    print(f"The desired maximum relative bias is {relative_bias_threshold*100}%.")
    print(f"This implies E[F] > 1 / {relative_bias_threshold}, so E[F] > {required_population_F}.\n")

    # --- Step 3: Relate E[F] to the non-centrality parameter (NCP or lambda) ---
    # E[F] = 1 + NCP / L, where L is the number of instruments.
    # We need to find the NCP that corresponds to the boundary condition E[F] = 10.
    # 10 = 1 + NCP / 1  => NCP = 9
    ncp_boundary = (required_population_F - 1) * num_instruments
    
    print("Step 2: Determine the corresponding non-centrality parameter (NCP)")
    print(f"The population F-statistic is related to the NCP by: E[F] = 1 + NCP / L.")
    print(f"For our boundary condition, the equation is: {required_population_F} = 1 + NCP / {num_instruments}.")
    print(f"Solving for the NCP gives a boundary value of {ncp_boundary}.\n")

    # --- Step 4: Find the critical value for the hypothesis test ---
    # We want to be 95% confident that the true NCP is greater than 9.
    # This is a hypothesis test: H0: NCP <= 9 vs. H1: NCP > 9.
    # We find the critical value that rejects H0 at a 5% significance level.
    # This is the 95th percentile of the F-statistic's distribution when NCP = 9.
    # For L=1, the F-statistic follows a non-central chi-squared distribution.
    df = num_instruments
    critical_value = scipy.stats.ncx2.ppf(confidence_level, df=df, nc=ncp_boundary)
    
    print("Step 3: Find the minimum sample F-statistic")
    print(f"To be {confidence_level*100}% confident that the true NCP is greater than {ncp_boundary},")
    print(f"we find the {int(confidence_level*100)}th percentile of the F-statistic's distribution assuming NCP = {ncp_boundary}.")
    print(f"For {num_instruments} instrument(s), this is a non-central chi-squared distribution with df={df} and ncp={ncp_boundary}.")
    print(f"The calculated critical value (minimum F-statistic) is: {critical_value:.4f}\n")

    # --- Step 5: Round to the nearest whole number ---
    final_answer = round(critical_value)
    
    print("Final Answer:")
    print(f"The minimum F-statistic required, rounded to the nearest whole number, is: {final_answer}")

if __name__ == '__main__':
    calculate_minimum_f_statistic()