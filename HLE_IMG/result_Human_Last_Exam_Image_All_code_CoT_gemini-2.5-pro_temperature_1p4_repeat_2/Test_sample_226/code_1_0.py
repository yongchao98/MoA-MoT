import numpy as np

def analyze_cell_counts():
    """
    This function analyzes and presents the cell count data for APT1 immunopositive cells
    based on the most plausible statement derived from the provided image.
    """
    # Data from statement A, which aligns with the visual evidence.
    # The values represent mean ± standard deviation for cells per mm^2.
    cell_counts = {
        "control": {"mean": 679.6, "std_dev": 59.32},
        "PD":      {"mean": 302.1, "std_dev": 111.5},
        "PDD":     {"mean": 283.2, "std_dev": 42.26}
    }

    # Print the findings in a clear format
    print("Analysis of APT1 Immunopositive Cell Counts:")
    print("="*45)
    print("The visual data suggests a significant decrease in APT1-positive cells in Parkinson's Disease (PD) and Parkinson's Disease with Dementia (PDD) brains compared to healthy controls.")
    print("The quantitative data that best represents this observation is:")

    # Iterate through the dictionary and print each group's data
    for group, data in cell_counts.items():
        mean = data['mean']
        std_dev = data['std_dev']
        print(f"- {group.upper()} group: {mean} ± {std_dev} cells per mm^2")

    print("\nThis indicates a reduction of over 50% in APT1-positive cell density in the disease states.")
    print("="*45)

# Execute the analysis
analyze_cell_counts()