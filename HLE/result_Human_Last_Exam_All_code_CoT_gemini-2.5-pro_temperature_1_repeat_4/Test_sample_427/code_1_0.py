import scipy.stats

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group for a two-sided
    Mann-Whitney U test to be significant after Bonferroni correction.
    """
    # --- Step 1: Define parameters ---
    initial_alpha = 0.05
    num_tests = 5
    bonferroni_alpha = initial_alpha / num_tests

    # The base patterns for the data in each group for n=5
    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [2, 2, 3, 3, 3]

    print(f"Goal: Find the minimum sample size 'n' per group to achieve a p-value < {bonferroni_alpha:.2f}\n")

    # --- Step 2: Iteratively find the minimum sample size ---
    # We will increase the sample size n in multiples of 5 to maintain the data pattern.
    # Start with a multiplier k=1 (n=5) and increase it.
    k = 1
    while True:
        n_per_group = 5 * k
        
        # Generate the data for the current sample size n by replicating the base pattern
        group1 = base_group1 * k
        group2 = base_group2 * k
        
        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            group1, group2, alternative='two-sided'
        )
        
        # --- Step 3: Check if the result is statistically significant ---
        if p_value < bonferroni_alpha:
            # We found the minimum n, so we print the final results and stop.
            print(f"Found: Minimum number of observations per group is {n_per_group}.")
            print("-" * 50)
            print("Final Statistically Significant Equation Details:")
            print(f"Bonferroni-corrected alpha: {initial_alpha} / {num_tests} = {bonferroni_alpha}")
            # The prompt asks to output each number in the final equation.
            # Here are the lists of numbers (observations) for each group.
            print(f"Group 1 (n={n_per_group}): {group1}")
            print(f"Group 2 (n={n_per_group}): {group2}")
            print(f"Result: Mann-Whitney U = {u_statistic:.1f}, p-value = {p_value:.5f}")
            print(f"The p-value ({p_value:.5f}) is less than the corrected alpha ({bonferroni_alpha}), achieving statistical significance.")
            
            # Return the final answer
            return n_per_group
            
        # If not significant, increase the multiplier for the next iteration
        k += 1

if __name__ == '__main__':
    find_minimum_sample_size()
