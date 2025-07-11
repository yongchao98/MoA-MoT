def analyze_flow_cytometry():
    """
    Analyzes the provided flow cytometry histograms to determine which antibody stained all cells.
    """
    print("Analyzing the flow cytometry histograms...")
    print("-" * 30)

    # Principle of interpretation
    print("Principle: If an antibody stains all cells, the entire cell population (red curve) will show increased fluorescence compared to the unstained control (black curve). This is seen as a complete shift of the population peak to the right.\n")

    # Analysis of Histogram A
    print("Analysis of Histogram A:")
    print("The red curve ('Antibody-stained cells') is a single peak, fully shifted to the right of the black curve ('Unstained cells').")
    print("Conclusion: Antibody A stained ALL cells in the sample.\n")
    stained_all_A = True

    # Analysis of Histogram B
    print("Analysis of Histogram B:")
    print("Similar to A, the red curve is a single peak completely shifted to a higher fluorescence.")
    print("Conclusion: Antibody B stained ALL cells in the sample.\n")
    stained_all_B = True

    # Analysis of Histogram C
    print("Analysis of Histogram C:")
    print("The red curve shows two distinct populations. One peak overlaps with the unstained control, and another is shifted to the right.")
    print("This indicates that some cells remained unstained.")
    print("Conclusion: Antibody C stained only a SUBSET of the cells, not all of them.\n")
    stained_all_C = False

    # Final Conclusion
    print("-" * 30)
    print("Final Conclusion:")
    print("The antibodies that stained all cells are those where the entire stained population shifted.")
    print(f"Antibody A stained all cells: {stained_all_A}")
    print(f"Antibody B stained all cells: {stained_all_B}")
    print(f"Antibody C stained all cells: {stained_all_C}")
    print("\nTherefore, antibodies A and B stained all cells in the sample.")

analyze_flow_cytometry()