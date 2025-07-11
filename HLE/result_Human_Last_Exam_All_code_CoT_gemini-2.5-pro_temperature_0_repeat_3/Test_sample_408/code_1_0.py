import numpy as np
from scipy import stats

def compare_extinction_points():
    """
    This function demonstrates the statistical procedure to compare the extinction points
    of two microbial cell types.
    
    1. It simulates extinction time data for two cell types.
    2. It checks if the data from both groups is normally distributed using the Shapiro-Wilk test.
    3. Based on the normality check, it selects the appropriate statistical test:
       - Unpaired t-test if data is normal.
       - Wilcoxon rank-sum (Mann-Whitney U) test if data is not normal.
    4. It performs the test and prints a conclusion.
    """
    # --- Step 1: Simulate Extinction Time Data ---
    # In a real scenario, you would use your own experimental data here.
    # We will simulate data for two cell types. Let's assume times are in hours.
    # We'll make the data non-normal to demonstrate the robust workflow.
    np.random.seed(42) # for reproducible results
    cell_type_A_extinction_times = stats.expon.rvs(loc=10, scale=5, size=15) # Non-normal data
    cell_type_B_extinction_times = stats.expon.rvs(loc=15, scale=5, size=17) # Non-normal data with a different mean

    print("Data for Cell Type A (extinction time in hours):")
    print([round(t, 2) for t in cell_type_A_extinction_times])
    print("\nData for Cell Type B (extinction time in hours):")
    print([round(t, 2) for t in cell_type_B_extinction_times])

    # Set the significance level (alpha)
    alpha = 0.05

    # --- Step 2: Check for Normality ---
    print("\n--- Checking for Normality (Shapiro-Wilk Test) ---")
    shapiro_A = stats.shapiro(cell_type_A_extinction_times)
    shapiro_B = stats.shapiro(cell_type_B_extinction_times)

    print(f"Cell Type A: Statistic={shapiro_A.statistic:.4f}, p-value={shapiro_A.pvalue:.4f}")
    print(f"Cell Type B: Statistic={shapiro_B.statistic:.4f}, p-value={shapiro_B.pvalue:.4f}")

    # Assume normality if p-value > alpha for both groups
    is_normal = shapiro_A.pvalue > alpha and shapiro_B.pvalue > alpha

    if is_normal:
        print("\nConclusion: Both datasets appear to be normally distributed.")
    else:
        print("\nConclusion: At least one dataset does not appear to be normally distributed.")

    # --- Step 3: Choose and Perform the Significance Test ---
    print("\n--- Performing Significance Test ---")
    if is_normal:
        print("Since data is normal, performing an Unpaired t-test...")
        test_statistic, p_value = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times, equal_var=True)
        print(f"T-test statistic: {test_statistic:.4f}")
        print(f"P-value: {p_value:.4f}")
    else:
        print("Since data is not normal, performing a Wilcoxon rank-sum test...")
        # Use alternative='two-sided' for a two-tailed test
        test_statistic, p_value = stats.mannwhitneyu(cell_type_A_extinction_times, cell_type_B_extinction_times, alternative='two-sided')
        print(f"Mann-Whitney U statistic: {test_statistic:.4f}")
        print(f"P-value: {p_value:.4f}")

    # --- Step 4: Interpret the Result ---
    print("\n--- Final Conclusion ---")
    if p_value < alpha:
        print(f"The p-value ({p_value:.4f}) is less than the significance level ({alpha}).")
        print("Therefore, we reject the null hypothesis.")
        print("There is a statistically significant difference between the extinction points of the two cell types.")
    else:
        print(f"The p-value ({p_value:.4f}) is greater than the significance level ({alpha}).")
        print("Therefore, we fail to reject the null hypothesis.")
        print("There is not enough evidence to claim a significant difference between the extinction points.")

if __name__ == '__main__':
    compare_extinction_points()
<<<D>>>