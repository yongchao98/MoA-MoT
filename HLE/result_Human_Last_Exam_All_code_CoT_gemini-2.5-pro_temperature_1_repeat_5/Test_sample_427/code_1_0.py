import sys
# It's good practice to install and import necessary libraries.
# If you don't have scipy, you can install it by running:
# pip install scipy
try:
    from scipy.stats import mannwhitneyu
except ImportError:
    print("Scipy library not found. Please install it using: pip install scipy")
    sys.exit(1)

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group for a two-sided
    Mann-Whitney U test to achieve significance after Bonferroni correction.
    """
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [3, 3, 3, 2, 2]

    # We start with one replication of the base data.
    replication_factor = 0
    
    # We will loop, increasing the sample size, until the p-value is significant.
    # A safety break is included to prevent an infinite loop.
    while replication_factor < 500:
        replication_factor += 1
        
        # The number of observations per group 'n' is the base size * replications
        n = len(base_group1) * replication_factor
        
        # Create the samples for this iteration by replicating the base data
        current_group1 = base_group1 * replication_factor
        current_group2 = base_group2 * replication_factor
        
        # Perform the two-sided Mann-Whitney U test
        # We use the default continuity correction which is standard.
        u_statistic, p_value = mannwhitneyu(current_group1, current_group2, alternative='two-sided')

        # Check if the p-value has reached the corrected significance level
        if p_value < corrected_alpha:
            print(f"To achieve statistical significance, we need a p-value < {corrected_alpha}")
            print("This is calculated by taking the original alpha and dividing by the number of tests.")
            print(f"Bonferroni Corrected Alpha = {original_alpha} / {num_tests} = {corrected_alpha}\n")
            print(f"The minimum number of observations per group required is: {n}")
            print(f"With n = {n}, the test yields a p-value of {p_value:.5f}.")
            print("The final comparison is:")
            print(f"{p_value:.5f} < {corrected_alpha}")
            return n

    print("Could not find a solution within the tested limits.")
    return None

if __name__ == '__main__':
    find_minimum_sample_size()
