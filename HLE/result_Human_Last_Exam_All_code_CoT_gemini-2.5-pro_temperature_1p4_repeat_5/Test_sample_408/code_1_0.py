import numpy as np
from scipy import stats
import math

def analyze_extinction_difference():
    """
    This function implements the correct statistical procedure to determine if the
    extinction points of two microbial cell types are significantly different.
    It follows these steps:
    1. Defines two sample datasets of extinction times (in hours).
    2. Checks each dataset for normality using the Shapiro-Wilk test.
    3. Based on the normality check, performs either:
       - An unpaired t-test (if data is normal).
       - A Wilcoxon rank-sum test (if data is not normal).
    4. Prints the details and a final conclusion.
    """
    # Step 1: Define the data.
    # This data represents the measured time-to-extinction from multiple replicates.
    cell_type_A_extinction_times = np.array([9.8, 10.5, 12.1, 8.4, 10.1, 9.3, 11.6, 7.9, 11.1, 10.5])
    cell_type_B_extinction_times = np.array([12.2, 13.6, 11.9, 14.1, 12.8, 13.2, 11.6, 13.0, 13.4, 12.5])

    print("This script determines if the extinction points of two cell types are significantly different.")
    print("\nStep 1: Define the data (time-to-extinction from replicates).")
    print(f"Cell Type A Extinction Times (hours): {cell_type_A_extinction_times}")
    print(f"Cell Type B Extinction Times (hours): {cell_type_B_extinction_times}\n")

    # Step 2: Check for normality using the Shapiro-Wilk test.
    # A p-value > 0.05 suggests the data can be treated as normally distributed.
    alpha = 0.05
    shapiro_stat_A, shapiro_p_A = stats.shapiro(cell_type_A_extinction_times)
    shapiro_stat_B, shapiro_p_B = stats.shapiro(cell_type_B_extinction_times)

    print("Step 2: Check if the data is normally distributed.")
    print(f"Cell A - Shapiro-Wilk Test: p-value = {shapiro_p_A:.4f}")
    print(f"Cell B - Shapiro-Wilk Test: p-value = {shapiro_p_B:.4f}")

    is_normal_A = shapiro_p_A > alpha
    is_normal_B = shapiro_p_B > alpha

    # Step 3: Perform the appropriate significance test.
    if is_normal_A and is_normal_B:
        print(f"\nResult: Both datasets appear normal (p > {alpha}).")
        print("Step 3: Perform an unpaired t-test to compare the means of the two groups.")

        t_statistic, p_value = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times, equal_var=False)

        # Outputting each number for the "final equation"
        mean_A = np.mean(cell_type_A_extinction_times)
        mean_B = np.mean(cell_type_B_extinction_times)
        var_A = np.var(cell_type_A_extinction_times, ddof=1) # ddof=1 for sample variance
        var_B = np.var(cell_type_B_extinction_times, ddof=1)
        n_A = len(cell_type_A_extinction_times)
        n_B = len(cell_type_B_extinction_times)
        
        # Manually calculate the denominator to show the steps
        denominator = math.sqrt(var_A/n_A + var_B/n_B)

        print("\n--- T-test 'Equation' Details ---")
        print("The t-statistic is calculated as: t = (mean1 - mean2) / sqrt(var1/n1 + var2/n2)")
        print(f"Mean of Cell A: {mean_A:.3f}")
        print(f"Mean of Cell B: {mean_B:.3f}")
        print(f"Variance of Cell A: {var_A:.3f}")
        print(f"Variance of Cell B: {var_B:.3f}")
        print(f"Sample size of Cell A: {n_A}")
        print(f"Sample size of Cell B: {n_B}")
        print(f"Calculated t-statistic = ({mean_A:.3f} - {mean_B:.3f}) / sqrt({var_A:.3f}/{n_A} + {var_B:.3f}/{n_B}) = {t_statistic:.4f}")
        print(f"P-value from test = {p_value:.6f}")

    else:
        print(f"\nResult: At least one dataset does not appear normal (p <= {alpha}).")
        print("Step 3: Perform a Wilcoxon rank-sum test (Mann-Whitney U).")
        
        u_statistic, p_value = stats.mannwhitneyu(cell_type_A_extinction_times, cell_type_B_extinction_times, alternative='two-sided')

        n_A = len(cell_type_A_extinction_times)
        n_B = len(cell_type_B_extinction_times)
        
        print("\n--- Mann-Whitney U Test 'Equation' Details ---")
        print("This test compares the sum of ranks for each group.")
        all_data = np.concatenate([cell_type_A_extinction_times, cell_type_B_extinction_times])
        ranks = stats.rankdata(all_data)
        rank_sum_A = np.sum(ranks[:n_A])
        u_A = rank_sum_A - (n_A * (n_A + 1)) / 2

        print(f"Sample size of Cell A (n1): {n_A}")
        print(f"Sample size of Cell B (n2): {n_B}")
        print(f"Sum of ranks for Cell A (R1): {rank_sum_A:.1f}")
        print(f"U statistic based on Cell A = {rank_sum_A:.1f} - ({n_A} * ({n_A} + 1) / 2) = {u_A:.1f}")
        print(f"Test Statistic from library (U): {u_statistic:.1f}")
        print(f"P-value from test = {p_value:.6f}")


    # Step 4: Draw a conclusion based on the p-value.
    print("\nStep 4: Draw a conclusion based on the p-value.")
    if p_value < alpha:
        print(f"Conclusion: Since the p-value ({p_value:.6f}) is less than {alpha}, we reject the null hypothesis.")
        print("The extinction points of the two cell types are significantly different.")
    else:
        print(f"Conclusion: Since the p-value ({p_value:.6f}) is not less than {alpha}, we fail to reject the null hypothesis.")
        print("There is not enough evidence to say the extinction points are significantly different.")

if __name__ == '__main__':
    analyze_extinction_difference()