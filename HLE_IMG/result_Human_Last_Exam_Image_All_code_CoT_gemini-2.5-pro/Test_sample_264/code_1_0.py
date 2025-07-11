def analyze_flow_cytometry():
    """
    This function analyzes the provided flow cytometry histograms to determine which antibody stained all cells.
    """
    print("--- Analysis of Flow Cytometry Data ---")
    print("The goal is to find the antibody that stained all cells in the sample.")
    print("This is indicated by a complete shift of the stained population (red curve) to a higher fluorescence compared to the unstained population (black curve), with no overlap.\n")

    # Analysis of Histogram A
    print("Analysis of A:")
    print("The red curve overlaps with the black curve. This means only a subpopulation of cells was stained.")
    print("Result for A: Did not stain all cells.\n")

    # Analysis of Histogram B
    print("Analysis of B:")
    print("The red curve is completely shifted to the right of the black curve, with no overlap.")
    print("This indicates that all cells were stained.")
    print("Result for B: Stained all cells.\n")

    # Analysis of Histogram C
    print("Analysis of C:")
    print("The red curve overlaps with the black curve. This means only a subpopulation of cells was stained.")
    print("Result for C: Did not stain all cells.\n")

    # Final Conclusion
    print("--- Conclusion ---")
    print("Only Antibody B stained all the cells in the sample.")

analyze_flow_cytometry()