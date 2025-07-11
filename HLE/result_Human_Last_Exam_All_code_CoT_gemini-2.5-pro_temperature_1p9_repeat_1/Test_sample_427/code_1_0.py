import scipy.stats as stats
import sys

# Define the base observations for each group
base_group1 = [1, 2, 3, 4, 5]
base_group2 = [2, 2, 3, 3, 3]

# Define the parameters for the statistical test
alpha = 0.05
num_tests = 5
corrected_alpha = alpha / num_tests

# We start by replicating the base groups once, then twice, and so on.
replication_factor = 1

print(f"Searching for the minimum sample size to achieve a p-value < {corrected_alpha}")
print("=" * 70)

while True:
    # Calculate the current sample size
    n_per_group = len(base_group1) * replication_factor
    
    # Create the groups for the current sample size by replication
    current_group1 = base_group1 * replication_factor
    current_group2 = base_group2 * replication_factor
    
    # Perform the two-sided Mann-Whitney U test
    # Set use_continuity=False for more accurate p-values with ties, like in this case.
    res = stats.mannwhitneyu(current_group1, current_group2, alternative='two-sided', use_continuity=False)
    p_value = res.pvalue

    print(f"Testing with n = {n_per_group} per group... P-value: {p_value:.6f}")
    
    # Check if the p-value is below the corrected alpha
    if p_value < corrected_alpha:
        print("\n" + "=" * 70)
        print("Statistical significance achieved!")
        print(f"The minimum number of observations per group is {n_per_group}.")
        print("\nFinal Equation Details:")
        print(f"Group 1 (n={n_per_group}) was constructed by repeating {base_group1}, {replication_factor} times.")
        print(f"Group 2 (n={n_per_group}) was constructed by repeating {base_group2}, {replication_factor} times.")
        print(f"The Mann-Whitney U statistic is {res.statistic:.1f}.")
        # Final output showing the numbers in the comparison
        print(f"The resulting p-value ({p_value:.6f}) is less than the Bonferroni-corrected alpha ({corrected_alpha}).")
        # Forcible flush to ensure the answer is printed last.
        sys.stdout.flush() 
        print(f"\n<<<{n_per_group}>>>")
        break
        
    # If not significant, increase the replication factor and try again
    replication_factor += 1
