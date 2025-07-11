import matplotlib.pyplot as plt
import numpy as np

def analyze_statements():
    """
    Analyzes the statements based on visual data and presents the quantitative
    data from statement A for comparison.
    """
    
    # Data from Answer Choice A
    groups = ['Control', 'PD', 'PDD']
    means = {
        'Control': 679.6,
        'PD': 302.1,
        'PDD': 283.2
    }
    std_devs = {
        'Control': 59.32,
        'PD': 111.5,
        'PDD': 42.26
    }

    print("--- Quantitative Data from Statement A ---")
    print("This statement claims the following cell densities:")
    for group in groups:
        # The prompt requested outputting each number in the final equation.
        # Here we print the values for each group as mean ± std.
        print(f"Group: {group}")
        print(f"Equation: Cell Density = {means[group]} ± {std_devs[group]} cells/mm²")
        print(f"Mean = {means[group]}")
        print(f"Standard Deviation = {std_devs[group]}")
        print("-" * 15)

    print("\n--- Visual Analysis Conclusion ---")
    print("Visual inspection of the images shows a clear INCREASE in the number and intensity")
    print("of APT1-positive cells from the 'control' group to the 'PD' and 'PDD' groups.")
    print("Statement A, which reports a DECREASE, is therefore inconsistent with the images.")
    print("Statement D, which reports an INCREASE in PDD brains, is consistent with the visual evidence.")

# Run the analysis
analyze_statements()
