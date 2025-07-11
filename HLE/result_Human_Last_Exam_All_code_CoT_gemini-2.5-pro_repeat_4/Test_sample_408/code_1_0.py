import numpy as np
from scipy import stats

def analyze_extinction_points():
    """
    This function demonstrates the statistical procedure to compare extinction points
    between two cell types.
    """
    # Step 1: Sample data representing extinction times (in hours) for two cell types.
    # In a real scenario, this data would come from your experiments.
    # We will make the data non-normal to demonstrate the Wilcoxon path.
    np.random.seed(0) # for reproducibility
    cell_type_A_extinction_times = np.random.exponential(scale=10, size=15) + 20
    cell_type_B_extinction_times = np.random.exponential(scale=15, size=16) + 25

    print("--- Extinction Point Data ---")
    print(f"Cell Type A (n={len(cell_type_A_extinction_times)}): {[round(x, 1) for x in cell_type_A_extinction_times]}")
    print(f"Cell Type B (n={len(cell_type_B_extinction_times)}): {[round(x, 1) for x in cell_type_B_extinction_times]}")
    print("-" * 30)

    # Step 2: Check for normality using the Shapiro-Wilk test
    alpha = 0.05
    stat_A, p_A = stats.shapiro(cell_type_A_extinction_times)
    stat_B, p_B = stats.shapiro(cell_type_B_extinction_times)

    print("--- Normality Check ---")
    print(f"Shapiro-Wilk test for Cell Type A: p-value = {p_A:.4f}")
    print(f"Shapiro-Wilk test for Cell Type B: p-value = {p_B:.4f}")

    # Step 3: Select and perform the appropriate significance test
    print("\n--- Significance Test ---")
    
    # Check if both p-values are greater than alpha
    if p_A > alpha and p_B > alpha:
        print(f"Data appears normally distributed (p > {alpha}).")
        print("Performing unpaired t-test...")

        # Perform t-test
        t_statistic, p_value = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times)

        print("\n--- T-test Results ---")
        print(f"Comparing samples: A vs B")
        print(f"T-statistic = {t_statistic:.4f}")
        print(f"P-value = {p_value:.4f}")

    else:
        print(f"At least one dataset does not appear normally distributed (p <= {alpha}).")
        print("Performing Wilcoxon rank-sum test (Mann-Whitney U)...")

        # Perform Wilcoxon rank-sum test
        u_statistic, p_value = stats.mannwhitneyu(cell_type_A_extinction_times, cell_type_B_extinction_times, alternative='two-sided')
        
        print("\n--- Wilcoxon Test Results ---")
        print(f"Comparing samples: A vs B")
        print(f"U-statistic = {u_statistic:.4f}")
        print(f"P-value = {p_value:.4f}")
        
    # Final conclusion
    print("\n--- Conclusion ---")
    if p_value < alpha:
        print(f"The p-value ({p_value:.4f}) is less than {alpha}, so we reject the null hypothesis.")
        print("The extinction points of the two cell types are significantly different.")
    else:
        print(f"The p-value ({p_value:.4f}) is not less than {alpha}, so we fail to reject the null hypothesis.")
        print("There is not enough evidence to claim the extinction points are significantly different.")

if __name__ == '__main__':
    analyze_extinction_points()