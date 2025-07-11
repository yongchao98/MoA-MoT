import numpy as np
from scipy import stats

def compare_extinction_points():
    """
    This function demonstrates the correct statistical procedure to compare the
    extinction points of two microbial cell types.
    
    1. It defines two sample datasets representing extinction points.
    2. It checks for normality using the Shapiro-Wilk test.
    3. Based on the normality check, it performs either an unpaired t-test (if normal)
       or a Wilcoxon rank-sum test (if not normal).
    4. It prints the results and a clear conclusion.
    """
    
    # --- Step 1: Define the extinction point data ---
    # These are sample extinction times (e.g., in hours) for multiple replicates
    # of two different cell types.
    # In a real experiment, you would replace this with your own data.
    cell_type_A_extinction = [20.5, 22.1, 21.6, 19.8, 20.9, 23.0, 21.2]
    cell_type_B_extinction = [23.5, 24.1, 19.9, 23.8, 25.0, 24.5, 23.9]
    
    print("Data for Cell Type A:", cell_type_A_extinction)
    print("Data for Cell Type B:", cell_type_B_extinction)
    print("-" * 30)

    # Significance level
    alpha = 0.05

    # --- Step 2: Check for normality using the Shapiro-Wilk test ---
    shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_extinction)
    shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_extinction)

    print("Checking for normality (Shapiro-Wilk Test):")
    print(f"P-value for Cell Type A: {shapiro_A_p:.4f}")
    print(f"P-value for Cell Type B: {shapiro_B_p:.4f}")
    
    is_normal = shapiro_A_p > alpha and shapiro_B_p > alpha
    
    print("-" * 30)
    
    # --- Step 3 & 4: Choose the correct test and interpret the results ---
    if is_normal:
        print("Conclusion: Both datasets appear to be normally distributed.")
        print("Performing an Unpaired T-test.")
        
        # Perform unpaired t-test
        t_stat, p_value = stats.ttest_ind(cell_type_A_extinction, cell_type_B_extinction)
        
        print(f"T-statistic: {t_stat:.4f}")
        print(f"P-value: {p_value:.4f}")
        
    else:
        print("Conclusion: At least one dataset does not appear to be normally distributed.")
        print("Performing a Wilcoxon Rank-Sum test.")
        
        # Perform Wilcoxon rank-sum test
        rank_stat, p_value = stats.ranksums(cell_type_A_extinction, cell_type_B_extinction)
        
        print(f"Test statistic: {rank_stat:.4f}")
        print(f"P-value: {p_value:.4f}")

    print("-" * 30)
    
    # Final conclusion based on the p-value
    if p_value < alpha:
        print(f"Final Result: The p-value ({p_value:.4f}) is less than {alpha}, so we reject the null hypothesis.")
        print("The difference in extinction points between the two cell types is statistically significant.")
    else:
        print(f"Final Result: The p-value ({p_value:.4f}) is not less than {alpha}, so we fail to reject the null hypothesis.")
        print("There is not enough evidence to claim a significant difference in extinction points.")

if __name__ == '__main__':
    compare_extinction_points()