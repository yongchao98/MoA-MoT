import numpy as np
from scipy import stats

def compare_extinction_points(data1, data2, alpha=0.05):
    """
    Checks for normality and performs the appropriate statistical test to compare two independent samples.

    Args:
        data1 (list or np.array): Extinction times for cell type 1.
        data2 (list or np.array): Extinction times for cell type 2.
        alpha (float): The significance level for the statistical tests.
    """
    print(f"Data for Cell Type 1: {data1}")
    print(f"Data for Cell Type 2: {data2}\n")

    # Step 1: Check for normality using the Shapiro-Wilk test
    print("--- Step 1: Checking for Normality ---")
    shapiro_stat1, shapiro_p1 = stats.shapiro(data1)
    shapiro_stat2, shapiro_p2 = stats.shapiro(data2)

    print(f"Cell Type 1: Shapiro-Wilk Test P-value = {shapiro_p1:.4f}")
    print(f"Cell Type 2: Shapiro-Wilk Test P-value = {shapiro_p2:.4f}")

    # Interpretation of normality test
    is_normal1 = shapiro_p1 > alpha
    is_normal2 = shapiro_p2 > alpha

    if is_normal1 and is_normal2:
        print("\nConclusion: Both datasets appear to be normally distributed.")
        test_type = "t-test"
    else:
        print("\nConclusion: At least one dataset does not appear to be normally distributed.")
        test_type = "wilcoxon"
    
    # Step 2: Perform the appropriate significance test
    print(f"\n--- Step 2: Performing {test_type.capitalize()} ---")
    
    if test_type == "t-test":
        # Perform independent (unpaired) t-test
        t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=False) # Welch's t-test is often preferred
        print(f"Independent t-test Results: T-statistic = {t_stat:.4f}, P-value = {p_value:.4f}")
    else:
        # Perform Wilcoxon rank-sum test (Mann-Whitney U)
        rank_stat, p_value = stats.ranksums(data1, data2)
        print(f"Wilcoxon Rank-Sum Test Results: Statistic = {rank_stat:.4f}, P-value = {p_value:.4f}")

    # Step 3: Interpret the final result
    print("\n--- Step 3: Final Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("Conclusion: There is a statistically significant difference in the extinction points between the two cell types.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is greater than or equal to {alpha}, we fail to reject the null hypothesis.")
        print("Conclusion: There is no statistically significant difference in the extinction points between the two cell types.")

# --- Example Usage ---

# Case 1: Normally distributed data
print("======================================================")
print("           CASE 1: NORMALLY DISTRIBUTED DATA          ")
print("======================================================")
# Sample data: extinction times in hours for two cell types from 10 replicates each
cell_type_A_normal = [10.2, 11.5, 9.8, 12.1, 10.5, 9.5, 11.8, 10.9, 12.5, 11.2]
cell_type_B_normal = [13.5, 14.1, 12.9, 14.5, 13.2, 14.8, 12.5, 13.9, 14.2, 13.0]
compare_extinction_points(cell_type_A_normal, cell_type_B_normal)


# Case 2: Non-normally distributed data
print("\n\n======================================================")
print("         CASE 2: NON-NORMALLY DISTRIBUTED DATA        ")
print("======================================================")
# Sample data (e.g., from a skewed distribution)
cell_type_C_non_normal = [8.1, 8.5, 8.7, 9.2, 9.5, 9.6, 10.1, 10.5, 15.2, 16.5]
cell_type_D_non_normal = [10.8, 11.2, 11.5, 12.0, 12.1, 12.3, 12.8, 17.5, 18.1, 18.5]
compare_extinction_points(cell_type_C_non_normal, cell_type_D_non_normal)