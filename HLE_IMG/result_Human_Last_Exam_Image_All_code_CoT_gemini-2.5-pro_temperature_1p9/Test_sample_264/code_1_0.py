def analyze_histograms():
    """
    Analyzes the provided flow cytometry histograms to determine which antibody stained all cells.
    """
    print("Analyzing the flow cytometry results:")
    print("-" * 35)

    # Analysis of Histogram A
    print("Analysis for Antibody A:")
    print("The red curve (stained) shows a peak shifted to the right, indicating positive staining.")
    print("However, this creates two distinct peaks when compared to the black curve (unstained).")
    print("This pattern means that only a sub-population of cells was stained.")
    print("Conclusion: Antibody A did NOT stain all cells.\n")

    # Analysis of Histogram B
    print("Analysis for Antibody B:")
    print("The red curve (stained) shows a single, uniform peak that is completely shifted to the right compared to the black curve (unstained).")
    print("There is no population remaining at the low fluorescence level of the unstained control.")
    print("This indicates that the entire cell population was stained.")
    print("Conclusion: Antibody B stained ALL cells.\n")

    # Analysis of Histogram C
    print("Analysis for Antibody C:")
    print("The red curve (stained) shows a major peak shifted to the right.")
    print("However, there is a prominent 'shoulder' on the left side of the red curve that overlaps with the unstained population's fluorescence level.")
    print("This shows that while most cells were stained, a small population remained unstained or was very dimly stained.")
    print("Conclusion: Antibody C did NOT stain all cells.\n")

    print("-" * 35)
    print("Final Answer:")
    print("Only antibody B resulted in the staining of the entire cell population.")

# Run the analysis
analyze_histograms()