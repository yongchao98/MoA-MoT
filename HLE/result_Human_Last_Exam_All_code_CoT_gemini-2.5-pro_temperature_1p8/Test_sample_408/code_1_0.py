import numpy as np
from scipy import stats

def compare_extinction_points():
    """
    This function demonstrates how to statistically compare the extinction points
    of two microbial cell types (A and B).

    It follows these steps:
    1. Defines sample data for extinction times for two cell types.
    2. Checks if the data for each cell type is normally distributed using the Shapiro-Wilk test.
    3. Based on the normality check, it chooses the appropriate statistical test:
       - Unpaired t-test if both samples are normal.
       - Wilcoxon rank-sum test if at least one sample is not normal.
    4. Performs the chosen test and prints the conclusion.
    """
    # Step 1: Sample data representing the time to extinction (e.g., in hours)
    # for multiple replicates of two different cell types.
    # In this example, Cell Type B's data is skewed to demonstrate the non-parametric path.
    cell_type_A_times = [48.2, 51.5, 53.1, 49.0, 47.7, 55.3]
    cell_type_B_times = [45.1, 42.8, 46.5, 38.2, 41.9, 25.0] # The value 25.0 makes this non-normal

    print(f"Data for Cell Type A: {cell_type_A_times}")
    print(f"Data for Cell Type B: {cell_type_B_times}\n")

    # Set the significance level (alpha)
    alpha = 0.05

    # Step 2: Check for normality using the Shapiro-Wilk test
    shapiro_A = stats.shapiro(cell_type_A_times)
    shapiro_B = stats.shapiro(cell_type_B_times)

    print(f"Normality Check (Shapiro-Wilk Test):")
    print(f"  - Cell Type A: P-value = {shapiro_A.pvalue:.4f}")
    print(f"  - Cell Type B: P-value = {shapiro_B.pvalue:.4f}\n")

    # Step 3: Decide which test to use
    is_A_normal = shapiro_A.pvalue > alpha
    is_B_normal = shapiro_B.pvalue > alpha

    if is_A_normal and is_B_normal:
        # --- The data appears to be normally distributed ---
        print("Data for both cell types appear normal. Performing an unpaired t-test.")
        
        # Perform unpaired t-test
        test_statistic, p_value = stats.ttest_ind(cell_type_A_times, cell_type_B_times, equal_var=True)

        print("\n--- Unpaired T-test Results ---")
        print(f"Equation Inputs (means): Mean(A) vs Mean(B)")
        print(f"Mean Extinction Time A: {np.mean(cell_type_A_times):.2f}")
        print(f"Mean Extinction Time B: {np.mean(cell_type_B_times):.2f}")
        print(f"T-statistic: {test_statistic:.4f}")

    else:
        # --- The data does not appear to be normally distributed ---
        print("At least one data set does not appear normal. Performing Wilcoxon rank-sum test.")

        # Perform Wilcoxon rank-sum test (Mann-Whitney U)
        test_statistic, p_value = stats.mannwhitneyu(cell_type_A_times, cell_type_B_times, alternative='two-sided')

        print("\n--- Wilcoxon Rank-Sum Test Results ---")
        print(f"Equation Inputs (medians): Median(A) vs Median(B)")
        print(f"Median Extinction Time A: {np.median(cell_type_A_times):.2f}")
        print(f"Median Extinction Time B: {np.median(cell_type_B_times):.2f}")
        print(f"U-statistic: {test_statistic:.4f}")
        
    # Step 4: Interpret the results
    print(f"P-value: {p_value:.4f}")
    
    if p_value < alpha:
        print(f"\nConclusion: Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("The difference in extinction points between the two cell types is statistically significant.")
    else:
        print(f"\nConclusion: Since the p-value ({p_value:.4f}) is not less than {alpha}, we fail to reject the null hypothesis.")
        print("There is not enough evidence to claim a significant difference in extinction points.")

# Run the analysis
if __name__ == "__main__":
    compare_extinction_points()
<<<D>>>