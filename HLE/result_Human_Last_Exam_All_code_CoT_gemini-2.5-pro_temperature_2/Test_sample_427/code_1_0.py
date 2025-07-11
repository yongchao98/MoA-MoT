import math
from scipy.stats import mannwhitneyu

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance with a Bonferroni-corrected p-value.
    """
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests
    
    print("This script finds the minimum sample size per group to achieve significance.")
    print("-" * 30)

    # Explain the Bonferroni Correction
    print("The Bonferroni correction is applied to the original p-value threshold.")
    print(f"Original alpha / Number of tests = Corrected alpha")
    print(f"{original_alpha} / {num_tests} = {corrected_alpha}")
    print(f"The new p-value required for significance is < {corrected_alpha}")
    print("-" * 30)

    # Start with n=5 and increment by 5
    n = 5
    # Set a safety limit to prevent an infinite loop
    max_n = 1000  

    while True:
        if n > max_n:
            print(f"Search stopped. No sample size found up to n = {max_n} that achieves the desired significance.")
            break

        # Generate Group 1 based on fixed proportions (1/5 for each score 1-5)
        count_per_score_g1 = n // 5
        group1 = []
        for score in [1, 2, 3, 4, 5]:
            group1.extend([score] * count_per_score_g1)

        # Generate Group 2 based on fixed proportions (2/5 for score 2, 3/5 for score 3)
        count_score_2_g2 = (n * 2) // 5
        count_score_3_g2 = (n * 3) // 5
        group2 = [2] * count_score_2_g2 + [3] * count_score_3_g2

        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = mannwhitneyu(group1, group2, alternative='two-sided')

        # Check if the result is significant
        if p_value < corrected_alpha:
            print(f"SUCCESS: Minimum sample size found.\n")
            print(f"Minimum observations per group (n): {n}\n")
            print(f"For n = {n}, the data for the two groups is:")
            print(f"Group 1 (n={len(group1)}): {group1}")
            print(f"Group 2 (n={len(group2)}): {group2}\n")
            
            print("With this sample size, the Mann-Whitney U test yields:")
            print(f"p-value = {p_value:.5f}\n")
            
            print("The result is statistically significant because:")
            print(f"p-value ({p_value:.5f}) < corrected alpha ({corrected_alpha})")
            
            # The final answer in the requested format
            print(f"\n<<<{n}>>>")
            break

        # If not significant, increment n and continue the loop
        n += 5

if __name__ == '__main__':
    find_minimum_sample_size()