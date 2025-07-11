import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group required to achieve
    statistical significance in a two-sided Mann-Whitney U test after Bonferroni correction.
    """
    # Define the base observations obtained in one session
    group1_base = [1, 2, 3, 4, 5]
    group2_base = [3, 3, 3, 2, 2]
    
    # Define the parameters for the statistical test
    initial_alpha = 0.05
    num_tests = 5
    
    # Apply Bonferroni correction
    corrected_alpha = initial_alpha / num_tests
    
    # We will check sample sizes in increments of 5 (one session's worth of data)
    # The loop will check for a number of sessions from 1 up to 100.
    for num_sessions in range(1, 101):
        
        # Calculate the total number of observations per group for the current iteration
        n_per_group = num_sessions * len(group1_base)
        
        # Create the full datasets by repeating the base observations
        current_group1 = group1_base * num_sessions
        current_group2 = group2_base * num_sessions
        
        # Perform the two-sided Mann-Whitney U test
        # 'auto' method handles ties and sample sizes appropriately
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            current_group1, 
            current_group2, 
            alternative='two-sided'
        )
        
        # Check if the p-value meets the corrected significance level
        if p_value < corrected_alpha:
            print(f"The Bonferroni corrected p-value threshold is {initial_alpha} / {num_tests} = {corrected_alpha}")
            print(f"Minimum number of observations per group found: {n_per_group}")
            print(f"At this sample size, the Mann-Whitney U test yields a p-value of {p_value:.6f}.")
            print("The final equation demonstrating significance is:")
            print(f"{p_value:.6f} < {corrected_alpha}")
            return n_per_group

    # This part would be reached if no solution is found within the loop's limit
    print("A sufficient sample size was not found within the tested range.")
    return None

# Run the function to find and print the result
find_minimum_observations()