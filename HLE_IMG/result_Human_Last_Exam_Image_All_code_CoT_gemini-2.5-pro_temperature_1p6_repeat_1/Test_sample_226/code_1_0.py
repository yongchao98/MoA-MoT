def analyze_apt1_data():
    """
    Presents the quantitative data for APT1 immunopositive cells
    based on the provided experimental image and options.
    """
    # Data from the most plausible answer choice (A)
    control_mean = 679.6
    control_sd = 59.32
    pd_mean = 302.1
    pd_sd = 111.5
    pdd_mean = 283.2
    pdd_sd = 42.26

    print("Analysis of APT1 Immunopositive Cells in Human Hippocampal Tissue:")
    print("-" * 70)
    print("The visual evidence in the images suggests a decrease in APT1-positive cells in the brains of patients with Parkinson's Disease (PD) and Parkinson's Disease with Dementia (PDD) compared to healthy controls.")
    print("The statement that aligns with this observation provides the following quantification:\n")

    print(f"Control (WT) Group:")
    print(f"\t- Mean cell count: {control_mean} cells/mm²")
    print(f"\t- Standard Deviation: {control_sd}\n")

    print(f"Parkinson's Disease (PD) Group:")
    print(f"\t- Mean cell count: {pd_mean} cells/mm²")
    print(f"\t- Standard Deviation: {pd_sd}\n")

    print(f"Parkinson's Disease with Dementia (PDD) Group:")
    print(f"\t- Mean cell count: {pdd_mean} cells/mm²")
    print(f"\t- Standard Deviation: {pdd_sd}\n")

    print("-" * 70)
    print("Conclusion: The quantitative data shows a significant decrease in APT1 immunopositive cells in both PD and PDD brains compared to control brains. This matches the visual evidence.")

# Run the analysis
analyze_apt1_data()
