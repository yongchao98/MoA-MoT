def analyze_flow_cytometry():
    """
    Analyzes the provided flow cytometry histograms to determine which antibody stained all cells.
    """
    print("To determine which antibody stained all cells, we need to find the histogram where the entire population of stained cells (red curve) has a higher fluorescence than the unstained cells (black curve).")
    print("\nAnalysis of each antibody:")
    print(" - Antibody A: In plot A, the red curve shows that some cells have become more fluorescent, but a significant population remains that overlaps with the unstained black curve. This means not all cells were stained.")
    print(" - Antibody B: In plot B, the entire red curve is shifted to the right of the black curve. There are no cells left in the position of the unstained control. This indicates that all cells in the sample were stained by antibody B.")
    print(" - Antibody C: In plot C, similar to plot A, the red curve shows a population of stained cells and another population that overlaps with the unstained control. This means not all cells were stained.")
    print("\nConclusion: Only antibody B stained all cells in the sample.")

analyze_flow_cytometry()