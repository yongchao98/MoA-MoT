import numpy as np
from scipy import stats

def simulate_cytokine_release():
    """
    Simulates and analyzes cytokine release from two groups of CAR T cells:
    1. Manufactured without Interleukin-15 (-IL15)
    2. Manufactured with Interleukin-15 (+IL15)

    The simulation is based on the prediction that the +IL15 group will show
    enhanced cytokine release.
    """
    print("--- Simulating a CAR T Cell Cytokine Release Experiment ---")

    # Seed for reproducible results
    np.random.seed(42)

    # Simulation Parameters
    sample_size = 12
    # Group 1: Without IL-15 (-IL15). Lower mean cytokine release (e.g., IFN-gamma pg/mL).
    mean_no_il15 = 450
    std_dev_no_il15 = 110
    # Group 2: With IL-15 (+IL15). Higher mean cytokine release.
    mean_with_il15 = 750
    std_dev_with_il15 = 130

    # Generate mock data based on a normal distribution
    group_no_il15 = np.random.normal(loc=mean_no_il15, scale=std_dev_no_il15, size=sample_size)
    group_with_il15 = np.random.normal(loc=mean_with_il15, scale=std_dev_with_il15, size=sample_size)

    # Ensure no negative cytokine values, which are biologically impossible
    group_no_il15[group_no_il15 < 0] = 0
    group_with_il15[group_with_il15 < 0] = 0

    print(f"\nGenerated Cytokine Data (IFN-gamma in pg/mL) for {sample_size} samples per group:")
    print("Group (-IL15):")
    # Outputting each number as requested
    print([round(x, 2) for x in group_no_il15])
    print("\nGroup (+IL15):")
    # Outputting each number as requested
    print([round(x, 2) for x in group_with_il15])

    # Calculate statistics for each group
    avg_no_il15 = np.mean(group_no_il15)
    avg_with_il15 = np.mean(group_with_il15)

    print("\n--- Analysis ---")
    print(f"Average Cytokine Release (-IL15): {avg_no_il15:.2f} pg/mL")
    print(f"Average Cytokine Release (+IL15): {avg_with_il15:.2f} pg/mL")
    
    # Perform an independent t-test to compare the means of the two groups.
    # The "final equation" here is the t-test. We will show the inputs (the data above)
    # and the resulting outputs (t-statistic and p-value).
    t_statistic, p_value = stats.ttest_ind(group_no_il15, group_with_il15, equal_var=False)

    print("\n--- Statistical Test (Independent T-test) ---")
    print(f"The statistical test compares the two means ({avg_no_il15:.2f} vs {avg_with_il15:.2f}) to see if the difference is likely real or due to random chance.")
    print(f"T-statistic: {t_statistic:.4f}")
    print(f"P-value: {p_value:.6f}")

    print("\n--- Conclusion ---")
    if p_value < 0.05:
        print("The p-value is less than 0.05, which is statistically significant.")
        print("This result supports the prediction that CAR T cells manufactured with IL-15 exhibit enhanced cytokine release.")
    else:
        print("The p-value is greater than 0.05. The difference in this simulation is not statistically significant.")

if __name__ == '__main__':
    simulate_cytokine_release()