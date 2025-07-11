import numpy as np
from scipy import stats

def compare_extinction_points():
    """
    This function demonstrates the statistical procedure to determine if the 
    extinction points of two microbial cell types are significantly different.
    """
    # Step 1: Define the data.
    # Assume we have collected the extinction times (in hours) for several 
    # replicate experiments for two independent cell types.
    cell_type_A_extinction_times = [22, 23, 23, 24, 25, 25, 26, 28]
    cell_type_B_extinction_times = [27, 28, 28, 29, 29, 30, 31, 33]

    print("This script determines if the extinction points of two cell types are significantly different.")
    print("\n--- Data ---")
    print(f"Extinction Times for Cell Type A: {cell_type_A_extinction_times}")
    print(f"Extinction Times for Cell Type B: {cell_type_B_extinction_times}")

    # Step 2: Check for normality using the Shapiro-Wilk test.
    # The null hypothesis of the Shapiro-Wilk test is that the data is normally distributed.
    alpha = 0.05
    shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_extinction_times)
    shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_extinction_times)

    print("\n--- Step 1: Normality Check (Shapiro-Wilk Test) ---")
    print(f"Cell Type A: p-value = {shapiro_A_p:.4f}")
    print(f"Cell Type B: p-value = {shapiro_B_p:.4f}")

    is_A_normal = shapiro_A_p > alpha
    is_B_normal = shapiro_B_p > alpha

    if is_A_normal and is_B_normal:
        print("\nConclusion: Both data sets appear to be normally distributed (p > 0.05).")
        # Step 3a: Perform an unpaired t-test if data is normal.
        print("Proceeding with an Unpaired T-test.")
        
        # Calculate components for the t-test equation
        mean_A = np.mean(cell_type_A_extinction_times)
        mean_B = np.mean(cell_type_B_extinction_times)
        std_A = np.std(cell_type_A_extinction_times, ddof=1)
        std_B = np.std(cell_type_B_extinction_times, ddof=1)
        n_A = len(cell_type_A_extinction_times)
        n_B = len(cell_type_B_extinction_times)
        
        print("\n--- Step 2a: Unpaired T-test Details ---")
        print("Equation components for t = (mean1 - mean2) / sqrt(s1^2/n1 + s2^2/n2)")
        print(f"Mean of Cell Type A (mean1): {mean_A:.2f}")
        print(f"Mean of Cell Type B (mean2): {mean_B:.2f}")
        print(f"Standard Deviation of A (s1): {std_A:.2f}")
        print(f"Standard Deviation of B (s2): {std_B:.2f}")
        print(f"Sample Size of A (n1): {n_A}")
        print(f"Sample Size of B (n2): {n_B}")
        
        t_stat, p_value = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times, equal_var=False)
        print("\n--- T-test Final Result ---")
        print(f"T-statistic = {t_stat:.4f}")

    else:
        print("\nConclusion: At least one data set is not normally distributed (p <= 0.05).")
        # Step 3b: Perform a Wilcoxon rank-sum test if data is not normal.
        print("Proceeding with the Wilcoxon rank-sum (Mann-Whitney U) test.")
        
        # The Mann-Whitney U test is the non-parametric equivalent of the t-test.
        u_stat, p_value = stats.mannwhitneyu(cell_type_A_extinction_times, cell_type_B_extinction_times, alternative='two-sided')
        print("\n--- Step 2b: Wilcoxon Rank-Sum Test Result ---")
        print(f"U-statistic = {u_stat:.1f}")

    # Step 4: Print the final conclusion.
    print(f"Final P-value = {p_value:.4f}")
    print("\n--- Final Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("The difference in extinction points between the two cell types is statistically significant.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is not less than {alpha}, we fail to reject the null hypothesis.")
        print("There is no statistically significant difference in extinction points between the two cell types.")

# Run the analysis
compare_extinction_points()