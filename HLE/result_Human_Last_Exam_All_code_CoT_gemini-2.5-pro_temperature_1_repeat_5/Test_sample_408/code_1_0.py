import numpy as np
from scipy import stats

def analyze_extinction_data(cell_type_1_data, cell_type_2_data, alpha=0.05):
    """
    Determines if the extinction points of two cell types are significantly different.
    
    This function follows the correct statistical procedure:
    1. It checks if the data for each cell type is normally distributed using the Shapiro-Wilk test.
    2. Based on the normality check, it selects the appropriate significance test:
       - Unpaired t-test if both datasets are normal.
       - Wilcoxon rank-sum test if at least one dataset is not normal.
    3. It prints the results and a conclusion.
    """
    print("--- Analysis of Extinction Point Data ---")
    # The data represents the measured extinction time (e.g., in hours) for each replicate.
    print(f"Data for Cell Type 1: {np.round(cell_type_1_data, 2)}")
    print(f"Data for Cell Type 2: {np.round(cell_type_2_data, 2)}\n")

    # Step 1: Check for normality in both datasets
    shapiro_stat1, norm_p_val1 = stats.shapiro(cell_type_1_data)
    shapiro_stat2, norm_p_val2 = stats.shapiro(cell_type_2_data)

    print("Step 1: Checking for Normality (Shapiro-Wilk Test)")
    print(f"Cell Type 1 normality p-value: {norm_p_val1:.4f}")
    print(f"Cell Type 2 normality p-value: {norm_p_val2:.4f}")

    is_normal1 = norm_p_val1 > alpha
    is_normal2 = norm_p_val2 > alpha

    # Step 2 & 3: Choose and perform the correct significance test
    if is_normal1 and is_normal2:
        print(f"\nConclusion: Both datasets appear normally distributed (p > {alpha}).")
        print("Performing an Unpaired t-test.\n")
        
        # Perform unpaired t-test
        t_statistic, p_value = stats.ttest_ind(cell_type_1_data, cell_type_2_data, equal_var=False) # Use Welch's t-test
        
        print("--- Unpaired t-test Results ---")
        print(f"T-statistic: {t_statistic:.4f}")
        print(f"P-value: {p_value:.4f}")

    else:
        print(f"\nConclusion: At least one dataset is not normally distributed (p <= {alpha}).")
        print("Performing a Wilcoxon rank-sum test.\n")

        # Perform Wilcoxon rank-sum (Mann-Whitney U) test
        u_statistic, p_value = stats.mannwhitneyu(cell_type_1_data, cell_type_2_data, alternative='two-sided')
        
        print("--- Wilcoxon Rank-Sum Test Results ---")
        print(f"U-statistic: {u_statistic:.4f}")
        print(f"P-value: {p_value:.4f}")

    # Step 4: Final conclusion based on the p-value
    print("\n--- Final Conclusion ---")
    if p_value < alpha:
        print(f"The result is significant (p < {alpha}). We conclude that the extinction points of the two cell types are significantly different.")
    else:
        print(f"The result is not significant (p >= {alpha}). We cannot conclude that the extinction points are different.")
    print("-" * 40 + "\n")


# --- Example 1: Using normally distributed data ---
print("Running Scenario 1: Data appears normal")
np.random.seed(42)
# Sample data: 12 replicates for each cell type. Times are in hours.
extinction_times_A = np.random.normal(loc=20, scale=2.5, size=12) # Mean 20h
extinction_times_B = np.random.normal(loc=24, scale=2.8, size=12) # Mean 24h
analyze_extinction_data(extinction_times_A, extinction_times_B)


# --- Example 2: Using non-normally distributed data ---
print("Running Scenario 2: Data does not appear normal")
np.random.seed(0)
# Sample data using a skewed distribution (log-normal)
extinction_times_C = np.random.lognormal(mean=3, sigma=0.2, size=15) # Median ~20h
extinction_times_D = np.random.lognormal(mean=3.4, sigma=0.25, size=15) # Median ~30h
analyze_extinction_data(extinction_times_C, extinction_times_D)
