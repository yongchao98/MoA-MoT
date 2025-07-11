import numpy as np
from scipy import stats

def analyze_extinction_times(data1, data2, group1_name="Cell Type A", group2_name="Cell Type B", alpha=0.05):
    """
    Performs a statistical analysis to compare extinction times of two cell types.
    
    1. Checks for normality in both data sets using the Shapiro-Wilk test.
    2. Based on normality, chooses and performs either an independent t-test (if normal)
       or a Mann-Whitney U test (if not normal).
    3. Prints the step-by-step analysis and conclusion.
    """
    print(f"--- Analysis of Extinction Times: {group1_name} vs {group2_name} ---")
    print(f"\nStep 1: Checking for Normality (using Shapiro-Wilk test at alpha={alpha})")
    
    # Perform Shapiro-Wilk test for normality on both datasets
    shapiro_stat1, shapiro_p1 = stats.shapiro(data1)
    shapiro_stat2, shapiro_p2 = stats.shapiro(data2)
    
    print(f"Normality test for {group1_name}: Statistic={shapiro_stat1:.4f}, p-value={shapiro_p1:.4f}")
    print(f"Normality test for {group2_name}: Statistic={shapiro_stat2:.4f}, p-value={shapiro_p2:.4f}")
    
    # The null hypothesis of the Shapiro-Wilk test is that the data is normal.
    # If p > alpha, we assume normality.
    is_normal1 = shapiro_p1 > alpha
    is_normal2 = shapiro_p2 > alpha
    
    if is_normal1 and is_normal2:
        print("\nConclusion: Both datasets appear to be normally distributed.")
        print("Step 2: Performing Independent t-test to compare the means.")
        
        # Calculate descriptive statistics (components of the t-test)
        mean1, mean2 = np.mean(data1), np.mean(data2)
        std1, std2 = np.std(data1, ddof=1), np.std(data2, ddof=1)
        n1, n2 = len(data1), len(data2)
        
        # Perform the t-test
        t_statistic, p_value = stats.ttest_ind(data1, data2, equal_var=False) # Welch's t-test is often preferred
        
        print("\n--- T-test Results ---")
        print(f"{group1_name}: Mean = {mean1:.2f}, SD = {std1:.2f}, n = {n1}")
        print(f"{group2_name}: Mean = {mean2:.2f}, SD = {std2:.2f}, n = {n2}")
        print(f"T-statistic = {t_statistic:.4f}")
        print(f"P-value = {p_value:.4f}")
        
    else:
        print("\nConclusion: At least one dataset does not appear to be normally distributed.")
        print("Step 2: Performing Mann-Whitney U test (non-parametric alternative).")
        
        # Perform the Mann-Whitney U test
        u_statistic, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        
        # Calculate medians for reporting
        median1, median2 = np.median(data1), np.median(data2)

        print("\n--- Mann-Whitney U Test Results ---")
        print(f"{group1_name}: Median = {median1:.2f}")
        print(f"{group2_name}: Median = {median2:.2f}")
        print(f"U-statistic = {u_statistic:.4f}")
        print(f"P-value = {p_value:.4f}")

    print("\nStep 3: Final Conclusion")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("There is a statistically significant difference in the extinction times between the two cell types.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is greater than or equal to {alpha}, we fail to reject the null hypothesis.")
        print("There is not enough evidence to claim a significant difference in the extinction times.")
    print("-" * 50 + "\n")


# --- Example 1: Data appears normal ---
print("Running Example 1 (Normal Data)")
# Sample data: Extinction times in hours for two cell types
cell_type_A_times = [20.5, 18.7, 21.2, 22.0, 19.5, 20.8, 23.1, 21.5]
cell_type_B_times = [24.1, 26.3, 25.5, 23.9, 27.0, 26.8, 24.5, 25.8]
analyze_extinction_times(cell_type_A_times, cell_type_B_times)


# --- Example 2: Data may not be normal ---
print("Running Example 2 (Non-Normal Data)")
# We add an outlier to one group to make it non-normal
cell_type_C_times = [15.1, 16.2, 14.8, 15.5, 16.0, 35.0] # The 35.0 is an outlier
cell_type_D_times = [18.2, 19.1, 19.5, 18.7, 20.1, 19.8]
analyze_extinction_times(cell_type_C_times, cell_type_D_times, "Cell Type C", "Cell Type D")