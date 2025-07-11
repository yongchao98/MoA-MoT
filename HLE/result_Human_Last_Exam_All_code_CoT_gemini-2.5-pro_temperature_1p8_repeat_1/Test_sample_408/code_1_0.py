import numpy as np
from scipy import stats

def compare_extinction_points(data1, data2, name1="Cell Type 1", name2="Cell Type 2", alpha=0.05):
    """
    Checks for normality and performs the appropriate statistical test
    to compare the extinction points of two cell types.

    Args:
    data1 (list or np.array): Extinction point data for the first cell type.
    data2 (list or np.array): Extinction point data for the second cell type.
    name1 (str): Name of the first cell type.
    name2 (str): Name of the second cell type.
    alpha (float): The significance level.
    """
    print(f"Comparing extinction points for '{name1}' and '{name2}'.\n")
    print(f"Data for {name1}: {data1}")
    print(f"Data for {name2}: {data2}")
    print("-" * 30)

    # Step 1: Check for normality using the Shapiro-Wilk test
    shapiro_stat1, shapiro_p1 = stats.shapiro(data1)
    shapiro_stat2, shapiro_p2 = stats.shapiro(data2)

    print("Step 1: Checking for normality with Shapiro-Wilk Test.")
    print(f"'{name1}': Statistics={shapiro_stat1:.3f}, p-value={shapiro_p1:.3f}")
    print(f"'{name2}': Statistics={shapiro_stat2:.3f}, p-value={shapiro_p2:.3f}")

    is_normal1 = shapiro_p1 > alpha
    is_normal2 = shapiro_p2 > alpha
    
    if is_normal1 and is_normal2:
        print("\nBoth datasets appear to be normally distributed (p > {alpha}).")
        # Step 2: If normal, perform an unpaired t-test
        print("Step 2: Performing an independent (unpaired) t-test.")
        ttest_stat, p_value = stats.ttest_ind(data1, data2)
        
        print("\n--- T-test Results ---")
        print(f"T-statistic: {ttest_stat:.3f}")

    else:
        print(f"\nAt least one dataset does not appear to be normally distributed (p <= {alpha}).")
        # Step 3: If not normal, perform the Wilcoxon rank-sum test (Mann-Whitney U)
        print("Step 2: Performing Wilcoxon rank-sum test (Mann-Whitney U).")
        # Note: 'alternative='two-sided' is the default
        mwu_stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')

        print("\n--- Wilcoxon Rank-Sum Test Results ---")
        print(f"U-statistic: {mwu_stat:.3f}")
        
    print(f"Final p-value: {p_value:.3f}")
    print("-" * 30)

    # Step 4: Interpret the result
    if p_value < alpha:
        print(f"Conclusion: The p-value is less than {alpha}, so we reject the null hypothesis.")
        print(f"The extinction points of '{name1}' and '{name2}' are significantly different.")
    else:
        print(f"Conclusion: The p-value is greater than or equal to {alpha}, so we fail to reject the null hypothesis.")
        print(f"There is no significant difference between the extinction points of '{name1}' and '{name2}'.")

# --- Example Usage ---

# Case 1: Both data sets are normally distributed
print("=========== CASE 1: NORMALLY DISTRIBUTED DATA ===========")
# Extinction time in hours for replicates of two cell types
cell_A_extinction_times = [48, 51, 52, 49, 47, 53, 50, 51]
cell_B_extinction_times = [55, 56, 54, 58, 57, 59, 56, 58]
compare_extinction_points(cell_A_extinction_times, cell_B_extinction_times, name1="Cell Type A", name2="Cell Type B")

print("\n\n=========== CASE 2: NON-NORMALLY DISTRIBUTED DATA ===========")
# Case 2: One data set is not normally distributed
cell_C_extinction_times = [40, 41, 42, 60, 62, 65, 66, 43] # Bimodal data
cell_D_extinction_times = [45, 46, 47, 49, 48, 45, 47, 46] # Normal data
compare_extinction_points(cell_C_extinction_times, cell_D_extinction_times, name1="Cell Type C", name2="Cell Type D")