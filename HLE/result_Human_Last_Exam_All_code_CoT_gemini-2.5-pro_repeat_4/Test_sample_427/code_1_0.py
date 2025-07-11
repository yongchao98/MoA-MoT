import scipy.stats

def solve_mann_whitney_power():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance for a Mann-Whitney U test after Bonferroni correction.
    """
    
    # --- Setup ---
    # The scores obtained in a single session for each group.
    group1_base = [1, 2, 3, 4, 5]
    group2_base = [3, 3, 3, 2, 2]
    
    # Significance level parameters
    original_alpha = 0.05
    num_tests = 5
    
    # --- Step 1: Calculate Bonferroni-corrected alpha ---
    corrected_alpha = original_alpha / num_tests
    
    # --- Step 2: Iteratively find the minimum sample size ---
    # 'k' represents the number of sessions run.
    k = 1
    while True:
        # Construct the full sample groups by repeating the base scores 'k' times.
        current_group1 = group1_base * k
        current_group2 = group2_base * k
        
        # Perform the two-sided Mann-Whitney U test.
        # Scipy automatically uses the normal approximation for larger sample sizes.
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            current_group1,
            current_group2,
            alternative='two-sided'
        )
        
        # Check if the p-value meets the corrected significance level.
        if p_value < corrected_alpha:
            n_per_group = len(current_group1)
            
            # --- Output the results ---
            print("To determine the minimum number of observations, we first adjust the p-value threshold using Bonferroni correction and then iteratively increase the sample size until we achieve a significant result.")
            print("\nStep 1: Calculate the corrected significance level (alpha)")
            print(f"Corrected Alpha = Original Alpha / Number of Tests")
            print(f"Corrected Alpha = {original_alpha} / {num_tests} = {corrected_alpha}\n")
            
            print(f"Step 2: Find the minimum number of sessions ('k') needed")
            print(f"By testing k=1, k=2, ..., we find that k={k} is the first point where the p-value is less than {corrected_alpha}.")
            print(f"Total observations per group (n) = 5 * k = {n_per_group}\n")
            
            print(f"Step 3: Present the final significant test result")
            print(f"Mann-Whitney U test with {n_per_group} observations per group yields a p-value of {p_value:.6f}.")
            print("The final equation comparing the p-value to the corrected alpha is:")
            print(f"{p_value:.6f} < {corrected_alpha}")
            
            print(f"\nTherefore, the minimum number of observations required per group is {n_per_group}.")
            
            # Final answer in the required format
            print(f"\n<<<{n_per_group}>>>")
            
            # Exit the loop once the answer is found.
            break
            
        # If not significant, increment the session count and try again.
        k += 1

# Run the solver function
solve_mann_whitney_power()