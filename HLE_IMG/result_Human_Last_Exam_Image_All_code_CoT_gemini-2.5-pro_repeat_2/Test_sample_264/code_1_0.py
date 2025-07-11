def analyze_flow_cytometry():
    """
    Analyzes the provided flow cytometry histograms to determine which antibody stained all cells.
    """
    print("Step-by-step Analysis of Flow Cytometry Histograms:")
    print("----------------------------------------------------\n")

    print("Understanding the Histograms:")
    print(" - The X-axis represents fluorescence intensity.")
    print(" - The Y-axis represents the count of cells.")
    print(" - The black curve is the 'unstained' control, showing baseline fluorescence (negative population).")
    print(" - The red curve is the 'antibody-stained' sample.")
    print(" - For an antibody to have stained ALL cells, the entire red curve must be shifted to the right (higher fluorescence) compared to the black curve, with no cells remaining in the negative/unstained region.\n")

    print("Analysis of each antibody:")
    print("1. Antibody A:")
    print("   The red curve shows a single population that is fully shifted to a higher fluorescence range compared to the unstained (black) control. This indicates all cells were stained.\n")

    print("2. Antibody B:")
    print("   The red curve is also completely shifted to the right, showing a clear separation from the unstained control. This indicates all cells were stained.\n")

    print("3. Antibody C:")
    print("   The red curve shows two populations: one population remains in the low fluorescence region (overlapping with the unstained control), and another population is shifted to the right. This means only a subset of cells was stained by antibody C.\n")

    print("Conclusion:")
    print("Antibodies A and B stained all cells in the sample, while antibody C stained only a fraction of the cells.")
    print("Therefore, the correct choice is A and B.")

analyze_flow_cytometry()