import numpy as np
from scipy import stats

def analyze_extinction_points():
    """
    This function demonstrates the statistical procedure for comparing the extinction points
    of two microbial cell types. It first checks for data normality and then selects the
    appropriate significance test (unpaired t-test for normal data, Wilcoxon rank-sum
    for non-normal data).
    """
    # --- Step 1: Define the data ---
    # These are hypothetical extinction times (e.g., in hours) from multiple
    # replicate experiments for two different cell types.
    # We will simulate a scenario where the data does not appear to be normal to show the full logic.
    print("--- Scenario: One dataset may not be normally distributed ---")

    cell_type_A_times = np.array([10.2, 11.1, 9.8, 10.5, 10.8, 11.5, 9.5])
    # Let's make Cell Type B's data skewed with a few high values (outliers)
    cell_type_B_times = np.array([11.8, 12.1, 11.5, 12.3, 15.0, 16.5, 11.9])

    print("Extinction Times for Cell Type A:", cell_type_A_times)
    print("Extinction Times for Cell Type B:", cell_type_B_times)
    print("-" * 30)

    # --- Step 2: Check for normality using the Shapiro-Wilk test ---
    # The null hypothesis is that the data is normally distributed.
    # We use a p-value threshold of 0.05.
    shapiro_A_stat, shapiro_A_pvalue = stats.shapiro(cell_type_A_times)
    shapiro_B_stat, shapiro_B_pvalue = stats.shapiro(cell_type_B_times)

    print(f"Shapiro-Wilk normality test for Cell A: p-value = {shapiro_A_pvalue:.4f}")
    print(f"Shapiro-Wilk normality test for Cell B: p-value = {shapiro_B_pvalue:.4f}")

    # --- Step 3: Choose and perform the appropriate significance test ---
    alpha = 0.05
    if shapiro_A_pvalue > alpha and shapiro_B_pvalue > alpha:
        # If both datasets appear normal (p > 0.05), use an unpaired t-test.
        print("\nBoth datasets appear to be normally distributed.")
        print("Performing an unpaired t-test...")
        
        t_stat, p_value = stats.ttest_ind(cell_type_A_times, cell_type_B_times)
        
        print(f"\nT-test Result: Statistic = {t_stat:.4f}, p-value = {p_value:.4f}")
        test_name = "t-test"

    else:
        # If at least one dataset does not appear normal (p <= 0.05), use a Wilcoxon rank-sum test.
        print("\nAt least one dataset does not appear to be normally distributed.")
        print("Performing a Wilcoxon rank-sum test...")

        ranksum_stat, p_value = stats.ranksums(cell_type_A_times, cell_type_B_times)
        
        print(f"\nWilcoxon Rank-Sum Test Result: Statistic = {ranksum_stat:.4f}, p-value = {p_value:.4f}")
        test_name = "Wilcoxon rank-sum test"

    # --- Step 4: Interpret the result ---
    print("\n--- Interpretation ---")
    if p_value < alpha:
        print(f"The p-value from the {test_name} is less than {alpha}.")
        print("Conclusion: The difference in extinction points between the two cell types is statistically significant.")
    else:
        print(f"The p-value from the {test_name} is not less than {alpha}.")
        print("Conclusion: There is not enough evidence to conclude a significant difference in extinction points.")

# Run the analysis
analyze_extinction_points()