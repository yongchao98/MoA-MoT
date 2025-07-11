import numpy as np
from scipy import stats

def analyze_extinction_times():
    """
    Analyzes the difference in extinction times between two cell types
    following the correct statistical procedure.
    """
    # Step 1: Define the data.
    # These are hypothetical extinction times (in hours) from 10 replicates for each cell type.
    # In a real scenario, you would replace this with your own measured data.
    cell_type_A_times = np.array([20.5, 21.0, 19.8, 20.2, 21.5, 22.0, 19.5, 20.8, 21.2, 20.0])
    cell_type_B_times = np.array([24.0, 24.5, 25.0, 23.8, 24.2, 26.0, 23.5, 25.5, 24.8, 25.2])

    print("--- Step 1: Data ---")
    print(f"Cell Type A Extinction Times (hours): {list(cell_type_A_times)}")
    print(f"Cell Type B Extinction Times (hours): {list(cell_type_B_times)}\n")

    # Step 2: Check for normality using the Shapiro-Wilk test.
    # The null hypothesis is that the data is normally distributed.
    # We use a significance level (alpha) of 0.05.
    alpha = 0.05
    shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_times)
    shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_times)

    print("--- Step 2: Normality Check (Shapiro-Wilk Test) ---")
    print(f"Cell Type A: p-value = {shapiro_A_p:.4f}")
    print(f"Cell Type B: p-value = {shapiro_B_p:.4f}")

    is_A_normal = shapiro_A_p > alpha
    is_B_normal = shapiro_B_p > alpha

    if is_A_normal and is_B_normal:
        print("Conclusion: Both samples appear to be normally distributed (p > 0.05).\n")
        test_to_use = 't-test'
    else:
        print("Conclusion: At least one sample does not appear to be normally distributed (p <= 0.05).\n")
        test_to_use = 'wilcoxon'

    # Step 3: Perform the appropriate significance test.
    print("--- Step 3: Significance Test ---")
    if test_to_use == 't-test':
        print("Data is normal. Performing Independent (Unpaired) T-test...")
        # Welch's t-test (equal_var=False) is a robust choice that doesn't assume equal variances.
        test_stat, p_value = stats.ttest_ind(cell_type_A_times, cell_type_B_times, equal_var=False)
        print(f"T-statistic = {test_stat:.4f}")
        print(f"P-value = {p_value:.4f}")
    else:
        print("Data is not normal. Performing Wilcoxon Rank-Sum (Mann-Whitney U) Test...")
        # The 'statistic' is the U statistic for the first sample (cell_type_A_times)
        test_stat, p_value = stats.mannwhitneyu(cell_type_A_times, cell_type_B_times, alternative='two-sided')
        print(f"U-statistic = {test_stat:.4f}")
        print(f"P-value = {p_value:.4f}")

    # Step 4: Interpret the final result.
    print("\n--- Step 4: Final Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("There IS a statistically significant difference between the extinction times of the two cell types.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is not less than {alpha}, we fail to reject the null hypothesis.")
        print("There is NOT a statistically significant difference between the extinction times of the two cell types.")

# Run the analysis
analyze_extinction_times()