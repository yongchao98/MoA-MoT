import numpy as np
from scipy import stats

def compare_extinction_points():
    """
    Compares the extinction points of two cell types to determine if they are
    significantly different, following a standard statistical procedure.
    """
    # --- Sample Data ---
    # These lists represent the measured extinction points (e.g., time in hours)
    # from multiple replicate experiments for two different cell types.
    # Replace this with your actual data.
    cell_type_A_extinction = [24.5, 25.1, 26.3, 23.9, 25.5, 24.8, 26.0, 25.2]
    cell_type_B_extinction = [27.1, 28.3, 27.9, 29.0, 28.5, 27.4, 28.1, 29.5]

    # Significance level (alpha)
    alpha = 0.05

    print("--- Step 1: Checking for Normality ---")
    # Perform the Shapiro-Wilk test for normality on each dataset.
    # The null hypothesis (H0) for this test is that the data is normally distributed.
    shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_extinction)
    shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_extinction)

    print(f"Normality check for Cell Type A: P-value = {shapiro_A_p:.4f}")
    print(f"Normality check for Cell Type B: P-value = {shapiro_B_p:.4f}")

    # --- Step 2: Select and Perform the Appropriate Significance Test ---
    # If p > alpha, we fail to reject H0, and assume the data is normal.
    if shapiro_A_p > alpha and shapiro_B_p > alpha:
        print("\nBoth datasets appear normally distributed. Performing Unpaired T-test.")
        
        # Perform the independent t-test
        t_stat, p_value = stats.ttest_ind(cell_type_A_extinction, cell_type_B_extinction, equal_var=True)
        
        print("\n--- Unpaired T-test Results ---")
        print(f"Data for Cell Type A: {cell_type_A_extinction}")
        print(f"Data for Cell Type B: {cell_type_B_extinction}")
        print(f"T-statistic = {t_stat:.4f}")
        print(f"P-value = {p_value:.4f}")

    else:
        print("\nAt least one dataset does not appear normal. Performing Wilcoxon Rank-Sum Test.")
        
        # Perform the Wilcoxon rank-sum test
        w_stat, p_value = stats.ranksums(cell_type_A_extinction, cell_type_B_extinction)

        print("\n--- Wilcoxon Rank-Sum Test Results ---")
        print(f"Data for Cell Type A: {cell_type_A_extinction}")
        print(f"Data for Cell Type B: {cell_type_B_extinction}")
        print(f"Wilcoxon statistic = {w_stat:.4f}")
        print(f"P-value = {p_value:.4f}")
        
    # --- Step 3: Interpret the Result ---
    print("\n--- Conclusion ---")
    if p_value < alpha:
        print(f"The p-value ({p_value:.4f}) is less than {alpha}.")
        print("We reject the null hypothesis: The difference in extinction points is statistically significant.")
    else:
        print(f"The p-value ({p_value:.4f}) is not less than {alpha}.")
        print("We fail to reject the null hypothesis: There is no significant statistical difference in extinction points.")

# Run the analysis
compare_extinction_points()