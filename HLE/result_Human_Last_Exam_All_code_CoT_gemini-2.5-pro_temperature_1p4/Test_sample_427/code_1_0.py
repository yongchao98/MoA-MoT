import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve statistical
    significance with Bonferroni correction for a Mann-Whitney U test.
    """
    original_alpha = 0.05
    num_tests = 5
    
    # Step 1: Calculate the Bonferroni-corrected p-value
    bonferroni_alpha = original_alpha / num_tests
    
    print(f"The original desired p-value was < {original_alpha}")
    print(f"The number of statistical tests is {num_tests}")
    print("To perform the Bonferroni correction, we divide the original p-value by the number of tests.")
    print("The new p-value threshold is:")
    print(f"{original_alpha} / {num_tests} = {bonferroni_alpha}")
    print("-" * 30)

    # Step 2: Define the base data for each group
    group1_base_scores = [1, 2, 3, 4, 5]
    group2_base_scores = [3, 3, 3, 2, 2]
    
    multiplier = 0
    p_value = 1.0 # Initialize p_value to be greater than alpha

    # Step 3 & 4: Iteratively increase sample size and test for significance
    while p_value >= bonferroni_alpha:
        multiplier += 1
        group1_sample = group1_base_scores * multiplier
        group2_sample = group2_base_scores * multiplier
        
        # Perform the two-sided Mann-Whitney U test
        # We use 'auto' for method which is recommended for larger samples with ties
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            group1_sample, 
            group2_sample, 
            alternative='two-sided',
            method='auto'
        )

    # Step 5: Output the result
    num_observations_per_group = len(group1_sample)
    
    print("Searching for the minimum number of observations per group...")
    print(f"Found it! When the number of observations per group is {num_observations_per_group}.")
    print(f"The resulting p-value is {p_value:.6f}, which is less than the corrected threshold of {bonferroni_alpha}.")
    print("\nTherefore, the minimum number of observations per group to achieve statistical significance is:")
    print(num_observations_per_group)

if __name__ == '__main__':
    find_minimum_observations()
    print("\n<<<45>>>")
