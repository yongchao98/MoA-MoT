import scipy.stats

def find_minimum_sample_size():
    """
    Calculates the minimum sample size per group for a two-sided Mann-Whitney U test
    to achieve significance after Bonferroni correction.
    """
    # Step 1: Define parameters and the corrected significance level
    original_alpha = 0.05
    num_tests = 5
    bonferroni_alpha = original_alpha / num_tests

    # Step 2: Define the base patterns for the two groups
    # Sorting group 2's pattern does not change the result but is good practice.
    pattern1 = [1, 2, 3, 4, 5]
    pattern2 = sorted([3, 3, 3, 2, 2])

    print(f"Original desired p-value: < {original_alpha}")
    print(f"Number of tests: {num_tests}")
    print(f"Bonferroni corrected p-value threshold: {original_alpha} / {num_tests} = {bonferroni_alpha}")
    print("\nSearching for the minimum number of observations (n) per group...")

    # Step 3: Iterate, starting from n=6 (since n=5 is given and not significant)
    n = 5
    while True:
        n += 1
        
        # Construct the groups for the current sample size 'n' by repeating the patterns
        # This approach ensures the distribution of values remains consistent as n grows.
        group1 = (pattern1 * (n // len(pattern1) + 1))[:n]
        group2 = (pattern2 * (n // len(pattern2) + 1))[:n]

        # Perform the two-sided Mann-Whitney U test
        statistic, p_value = scipy.stats.mannwhitneyu(group1, group2, alternative='two-sided')

        # Check if the result is statistically significant
        if p_value < bonferroni_alpha:
            print(f"\nFound significance at n = {n}.")
            print("The final comparison is:")
            # Step 4: Output the final equation and the result
            print(f"p-value ({p_value:.6f}) < Bonferroni-corrected alpha ({bonferroni_alpha})")
            
            print(f"\nThe minimum number of observations per group to achieve statistical significance is {n}.")
            return n

# Run the function to find the answer
minimum_n = find_minimum_sample_size()
print(f"<<<{minimum_n}>>>")