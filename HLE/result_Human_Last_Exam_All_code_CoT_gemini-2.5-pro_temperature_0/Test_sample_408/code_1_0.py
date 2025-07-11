import numpy as np
from scipy import stats

def analyze_extinction_points():
    """
    This function demonstrates the statistical procedure to determine if the 
    extinction points of two microbial cell types are significantly different.
    
    The logic is as follows:
    1. Define sample data for extinction times for two cell types.
    2. Check if the data in each group is normally distributed using the Shapiro-Wilk test.
    3. If both are normal, perform an independent t-test.
    4. If at least one is not normal, perform a Wilcoxon rank-sum (Mann-Whitney U) test.
    5. Print the steps, results, and a final conclusion.
    """
    # For reproducibility of the random data
    np.random.seed(42)

    # Step 1: Sample Data
    # Let's assume we have extinction times (in hours) from 12 replicates for each cell type.
    # We create two sets of data to compare.
    cell_type_A_times = np.random.normal(loc=48, scale=6, size=12)  # Mean extinction time ~48h
    cell_type_B_times = np.random.normal(loc=57, scale=7, size=12)  # Mean extinction time ~57h

    print("--- Comparing Extinction Points of Two Cell Types ---")
    print("\nExtinction Times for Cell Type A (hours):")
    # The numbers in the "equation" are the raw data points for each group
    print(' '.join([f'{x:.1f}' for x in cell_type_A_times]))
    
    print("\nExtinction Times for Cell Type B (hours):")
    print(' '.join([f'{x:.1f}' for x in cell_type_B_times]))
    print("-" * 55)

    # Set the significance level (alpha)
    alpha = 0.05

    # Step 2: Check for Normality
    shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_times)
    shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_times)

    print(f"\nStep 1: Checking for Normality (alpha = {alpha})")
    print(f"Shapiro-Wilk test for Cell Type A: p-value = {shapiro_A_p:.4f}")
    print(f"Shapiro-Wilk test for Cell Type B: p-value = {shapiro_B_p:.4f}")

    is_A_normal = shapiro_A_p > alpha
    is_B_normal = shapiro_B_p > alpha

    # Step 3 & 4: Select and Perform the appropriate test
    if is_A_normal and is_B_normal:
        print("\nConclusion: Both datasets appear normally distributed. Proceeding with t-test.")
        print("-" * 55)
        
        # Perform independent t-test
        t_statistic, p_value = stats.ttest_ind(cell_type_A_times, cell_type_B_times)
        
        print("Step 2: Independent t-test Results")
        # The numbers in the final "equation" are the test statistic and the p-value
        print(f"t-statistic = {t_statistic:.4f}")
        print(f"p-value = {p_value:.4f}")

    else:
        print("\nConclusion: At least one dataset is not normally distributed. Proceeding with Wilcoxon test.")
        print("-" * 55)
        
        # Perform Wilcoxon rank-sum test
        u_statistic, p_value = stats.mannwhitneyu(cell_type_A_times, cell_type_B_times, alternative='two-sided')
        
        print("Step 2: Wilcoxon Rank-Sum Test Results")
        # The numbers in the final "equation" are the test statistic and the p-value
        print(f"U-statistic = {u_statistic:.1f}")
        print(f"p-value = {p_value:.4f}")

    # Step 5: Interpret the final result
    print("-" * 55)
    print("Final Conclusion:")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("The difference in extinction times between the two cell types is statistically significant.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is not less than {alpha}, we fail to reject the null hypothesis.")
        print("The difference in extinction times between the two cell types is not statistically significant.")

# Run the analysis
analyze_extinction_points()