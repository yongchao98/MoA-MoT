import scipy.stats as stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve statistical
    significance with a Bonferroni-corrected Mann-Whitney U test.
    """
    # Step 1: Define initial parameters and calculate the Bonferroni-corrected alpha
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    print("--- Bonferroni Correction ---")
    # The final equation as requested
    print(f"Original alpha / Number of tests = Corrected alpha")
    print(f"{original_alpha} / {num_tests} = {corrected_alpha}")
    print(f"The new p-value threshold is {corrected_alpha}.\n")

    # Step 2: Define the base patterns for the two groups
    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [3, 3, 3, 2, 2]
    
    # Step 3: Iterate, increasing sample size until significance is found
    print("--- Searching for Minimum Sample Size ---")
    replication_factor = 1
    while True:
        n_per_group = len(base_group1) * replication_factor
        
        # Create the groups by replicating the base patterns
        current_group1 = base_group1 * replication_factor
        current_group2 = base_group2 * replication_factor
        
        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = stats.mannwhitneyu(
            current_group1, 
            current_group2, 
            alternative='two-sided'
        )
        
        print(f"Testing with n = {n_per_group} per group... p-value = {p_value:.6f}")
        
        # Check if the result is statistically significant
        if p_value < corrected_alpha:
            print("\n--- Result ---")
            print(f"Significance achieved at n = {n_per_group} observations per group.")
            print(f"The p-value ({p_value:.6f}) is less than the corrected alpha ({corrected_alpha}).")
            print(f"Therefore, the minimum number of observations per group required is {n_per_group}.")
            return n_per_group
            
        replication_factor += 1
        
        # Safety break to prevent an infinite loop
        if replication_factor > 100:
            print("Could not find a solution within a reasonable number of iterations.")
            return None

# Run the function to find the answer
minimum_n = find_minimum_observations()
print(f'<<<{minimum_n}>>>')
