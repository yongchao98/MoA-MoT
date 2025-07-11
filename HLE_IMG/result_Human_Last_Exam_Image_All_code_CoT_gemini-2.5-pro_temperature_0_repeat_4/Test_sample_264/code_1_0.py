def analyze_flow_cytometry():
    """
    Analyzes conceptual flow cytometry data to determine which antibody stained all cells.
    """
    print("Principle of Analysis:")
    print("In flow cytometry, a histogram shows cell count versus fluorescence intensity.")
    print("The black curve represents unstained cells (negative control), showing baseline fluorescence.")
    print("The red curve represents antibody-stained cells.")
    print("If an antibody stains ALL cells, the entire red curve should shift to a higher fluorescence range compared to the black curve, with no cells remaining at the baseline level.\n")

    print("--- Analysis of Each Antibody ---")

    # Analysis of Histogram A
    print("Histogram A:")
    print("The red curve shows two populations: one with high fluorescence (stained) and one that overlaps with the black curve (unstained).")
    print("Conclusion: Antibody A did NOT stain all cells.\n")

    # Analysis of Histogram B
    print("Histogram B:")
    print("The red curve is completely shifted to the right of the black curve. There is no population of cells left at the low, unstained fluorescence level.")
    print("Conclusion: Antibody B stained ALL cells in the sample.\n")

    # Analysis of Histogram C
    print("Histogram C:")
    print("The red curve shows a major stained population, but also a 'shoulder' on the left that overlaps with the black unstained curve.")
    print("Conclusion: Antibody C did NOT stain all cells.\n")

    # Final Conclusion
    print("--- Final Result ---")
    print("Only antibody B caused a uniform shift in fluorescence for the entire cell population.")
    print("Therefore, antibody B is the one that stained all cells in the sample.")

analyze_flow_cytometry()
<<<E>>>