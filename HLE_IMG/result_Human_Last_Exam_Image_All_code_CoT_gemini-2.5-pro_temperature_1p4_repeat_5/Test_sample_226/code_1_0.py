import sys

# Suppress any potential warning about running as a script
if "ipykernel" in sys.modules:
    pass

def analyze_cell_counts():
    """
    This function analyzes the cell count data from the most plausible statement
    and calculates the percentage change relative to the control group.
    """
    # Data from statement A (cells per mm2)
    control_count = 679.6
    pd_count = 302.1
    pdd_count = 283.2

    print("Analyzing the quantitative data from the most likely true statement (A)...")
    print("-" * 30)
    print(f"Control group (WT) mean cell count: {control_count}")
    print(f"Parkinson's Disease (PD) group mean cell count: {pd_count}")
    print(f"Parkinson's with Dementia (PDD) group mean cell count: {pdd_count}")
    print("-" * 30)

    # Calculate the percentage decrease compared to the control
    decrease_pd = ((control_count - pd_count) / control_count) * 100
    decrease_pdd = ((control_count - pdd_count) / control_count) * 100

    print("This corresponds to a significant reduction in APT1 immunopositive cells:")
    print(f"Decrease in PD brains compared to control: {decrease_pd:.2f}%")
    print(f"Decrease in PDD brains compared to control: {decrease_pdd:.2f}%")
    print("\nThis quantitative analysis confirms the visual observation that APT1 expression is reduced in PD and PDD brains.")

analyze_cell_counts()