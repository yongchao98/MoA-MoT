import scipy.stats as stats
import numpy as np

def analyze_extinction_points():
    """
    Analyzes the difference in extinction points between two cell types.

    This function follows these steps:
    1. Defines sample data for extinction times of two cell types.
    2. Checks if the data in each group is normally distributed using the Shapiro-Wilk test.
    3. Based on the normality check, it selects and performs the appropriate statistical test:
       - Unpaired t-test if data is normal.
       - Wilcoxon rank-sum (Mann-Whitney U) test if data is not normal.
    4. Prints the results and a clear conclusion.
    """
    # Sample data: Extinction times (e.g., in hours) for two cell types.
    # Each list represents multiple experimental replicates.
    cell_type_A_extinction_times = np.array([20.5, 22.1, 21.8, 23.0, 21.5, 22.5, 20.9, 22.8])
    cell_type_B_extinction_times = np.array([24.1, 25.5, 26.2, 23.9, 25.8, 26.5, 24.7, 25.1])

    print("--- Data ---")
    print(f"Cell Type A Extinction Times: {cell_type_A_extinction_times}")
    print(f"Cell Type B Extinction Times: {cell_type_B_extinction_times}")
    print("-" * 20)

    # Significance level
    alpha = 0.05

    # Step 1: Check for normality using the Shapiro-Wilk test
    print("\n--- Step 1: Normality Check ---")
    shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_extinction_times)
    shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_extinction_times)

    print(f"Shapiro-Wilk test for Cell Type A: Statistic={shapiro_A_stat:.4f}, p-value={shapiro_A_p:.4f}")
    print(f"Shapiro-Wilk test for Cell Type B: Statistic={shapiro_B_stat:.4f}, p-value={shapiro_B_p:.4f}")

    # Assume normality if p-value > alpha for both datasets
    is_normal = shapiro_A_p > alpha and shapiro_B_p > alpha

    # Step 2: Perform the appropriate significance test
    print("\n--- Step 2: Significance Test ---")
    if is_normal:
        print("Data from both groups appear to be normally distributed (p > 0.05).")
        print("Performing an unpaired t-test.")
        
        # Perform unpaired t-test
        t_stat, p_value = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times)
        
        print("\n--- Test Results ---")
        print(f"T-statistic = {t_stat:.4f}")
        print(f"P-value = {p_value:.4f}")

    else:
        print("Data from at least one group does not appear to be normally distributed (p <= 0.05).")
        print("Performing a Wilcoxon rank-sum test (Mann-Whitney U).")
        
        # Perform Wilcoxon rank-sum test
        u_stat, p_value = stats.mannwhitneyu(cell_type_A_extinction_times, cell_type_B_extinction_times, alternative='two-sided')
        
        print("\n--- Test Results ---")
        print(f"U-statistic = {u_stat:.4f}")
        print(f"P-value = {p_value:.4f}")

    # Step 3: Interpret the results
    print("\n--- Step 3: Conclusion ---")
    if p_value < alpha:
        print(f"The p-value ({p_value:.4f}) is less than our significance level ({alpha}).")
        print("Therefore, we reject the null hypothesis.")
        print("Conclusion: There is a statistically significant difference between the extinction points of the two cell types.")
    else:
        print(f"The p-value ({p_value:.4f}) is not less than our significance level ({alpha}).")
        print("Therefore, we fail to reject the null hypothesis.")
        print("Conclusion: There is no statistically significant difference between the extinction points of the two cell types.")

# Run the analysis
analyze_extinction_points()