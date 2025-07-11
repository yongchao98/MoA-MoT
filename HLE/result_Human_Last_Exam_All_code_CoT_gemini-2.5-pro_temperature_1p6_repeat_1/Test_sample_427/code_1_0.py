import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve statistical
    significance with a Bonferroni-corrected p-value in a Mann-Whitney U test.
    """
    # Step 1: Define initial parameters and calculate the corrected alpha
    pattern_group1 = [1, 2, 3, 4, 5]
    pattern_group2 = [3, 3, 3, 2, 2]
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    print("--- Bonferroni Correction Calculation ---")
    print(f"Original alpha: {original_alpha}")
    print(f"Number of tests: {num_tests}")
    # The final equation as requested by the prompt
    print(f"Corrected alpha = {original_alpha} / {num_tests} = {corrected_alpha}")
    print("-" * 40)

    # Step 2: Iterate to find the minimum sample size
    num_sessions = 0
    min_observations = 0

    while True:
        num_sessions += 1
        n_per_group = len(pattern_group1) * num_sessions

        # Create sample data by replicating the initial patterns
        group1_data = pattern_group1 * num_sessions
        group2_data = pattern_group2 * num_sessions

        # Perform the two-sided Mann-Whitney U test
        # method='auto' correctly handles ties and sample size
        res = scipy.stats.mannwhitneyu(group1_data, group2_data, alternative='two-sided')
        p_value = res.pvalue

        print(f"Testing with n = {n_per_group} observations per group...")
        print(f"Calculated p-value: {p_value:.6f}")

        # Step 3: Check if the p-value is below the corrected alpha
        if p_value < corrected_alpha:
            print(f"\nSuccess! The p-value ({p_value:.6f}) is less than the corrected alpha ({corrected_alpha}).")
            min_observations = n_per_group
            break
        else:
            print(f"The p-value is not less than {corrected_alpha}. Increasing sample size.\n")
    
    print("-" * 40)
    print(f"The minimum number of observations per group required to achieve statistical significance is {min_observations}.")


if __name__ == '__main__':
    find_minimum_observations()