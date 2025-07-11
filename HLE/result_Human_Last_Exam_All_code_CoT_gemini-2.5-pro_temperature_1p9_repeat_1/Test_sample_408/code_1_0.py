import numpy as np
from scipy import stats

def analyze_extinction_points():
    """
    Demonstrates the statistical procedure to compare extinction points
    of two microbial cell types.
    """
    # --- Step 0: Create Sample Data ---
    # For a real scenario, you would load your own data here.
    # Let's generate some sample data for demonstration. We will run two scenarios:
    # 1. Both datasets are normal, leading to a t-test.
    # 2. One dataset is not normal, leading to a Wilcoxon test.

    print("=========================================================")
    print("SCENARIO 1: Both distributions appear normally distributed")
    print("=========================================================\n")
    # Use a seed for reproducibility of the random data
    np.random.seed(42)
    # Extinction times (e.g., in hours) for 30 replicates of each cell type.
    cell_type_A_normal = np.random.normal(loc=10.0, scale=1.5, size=30)
    cell_type_B_normal = np.random.normal(loc=11.5, scale=1.7, size=30)
    
    # Run the analysis on the first scenario
    run_comparison(cell_type_A_normal, cell_type_B_normal)

    print("\n\n=========================================================")
    print("SCENARIO 2: One distribution is not normally distributed")
    print("=========================================================\n")
    np.random.seed(0)
    # One normal dataset, and one non-normal (exponential) dataset
    cell_type_A_normal = np.random.normal(loc=10.0, scale=1.5, size=30)
    cell_type_B_non_normal = np.random.exponential(scale=11.0, size=30)
    
    # Run the analysis on the second scenario
    run_comparison(cell_type_A_non_normal, cell_type_B_non_normal)


def run_comparison(data_A, data_B):
    """
    Performs the full analysis for two given datasets.
    """
    print("--- Step 1: Raw Data ---")
    print(f"Extinction times for Cell Type A: \n{np.round(data_A, 1)}")
    print(f"\nExtinction times for Cell Type B: \n{np.round(data_B, 1)}")

    # Set the significance level (alpha)
    alpha = 0.05
    print(f"\n--- Step 2: Normality Check (alpha = {alpha}) ---")

    # Perform Shapiro-Wilk test on both datasets
    shapiro_stat_A, shapiro_p_A = stats.shapiro(data_A)
    shapiro_stat_B, shapiro_p_B = stats.shapiro(data_B)

    print(f"Cell Type A: Shapiro-Wilk p-value = {shapiro_p_A:.4f}")
    print(f"Cell Type B: Shapiro-Wilk p-value = {shapiro_p_B:.4f}")

    is_A_normal = shapiro_p_A > alpha
    is_B_normal = shapiro_p_B > alpha

    # --- Step 3: Select and Perform Significance Test ---
    if is_A_normal and is_B_normal:
        print("\nDecision: Both datasets appear normal. Using Unpaired T-test.")
        test_statistic, p_value = stats.ttest_ind(data_A, data_B)
        
        print("\n--- T-test Results ---")
        # To "output each number in the final equation", we show the key components
        print(f"Mean of A: {np.mean(data_A):.2f}")
        print(f"Mean of B: {np.mean(data_B):.2f}")
        print(f"t-statistic: {test_statistic:.4f}")
        print(f"p-value: {p_value:.4f}")

    else:
        print("\nDecision: At least one dataset is not normal. Using Wilcoxon rank-sum test.")
        statistic, p_value = stats.ranksums(data_A, data_B)
        
        print("\n--- Wilcoxon Rank-Sum Test Results ---")
        print(f"Test statistic (z-score): {statistic:.4f}")
        print(f"p-value: {p_value:.4f}")

    # --- Step 4: Interpret the Result ---
    print("\n--- Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than alpha ({alpha}), we reject the null hypothesis.")
        print("This suggests a SIGNIFICANT difference between the extinction points of the two cell types.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is greater than alpha ({alpha}), we fail to reject the null hypothesis.")
        print("There is NOT enough evidence to claim a significant difference between the extinction points.")

# Run the full analysis demonstration
analyze_extinction_points()
