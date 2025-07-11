def analyze_flow_cytometry():
    """
    Analyzes the provided flow cytometry histograms to determine which antibody
    stained all cells in the sample.
    """
    print("Analyzing the flow cytometry results for antibodies A, B, and C.")
    print("The goal is to find the antibody that stained ALL cells.")
    print("This means the red curve (stained sample) should be completely shifted to the right of the black curve (unstained control), with no overlap.\n")

    # Analysis of Histogram A
    print("--- Analysis of Histogram A ---")
    print("The red curve is shifted to the right, indicating staining.")
    print("However, there is an overlap between the black and red curves.")
    print("This means some cells in the stained sample have fluorescence levels similar to unstained cells.")
    print("Conclusion: Antibody A did NOT stain all cells.\n")

    # Analysis of Histogram B
    print("--- Analysis of Histogram B ---")
    print("The red curve is completely separated from the black curve.")
    print("There is a clear gap between the two populations, with the red curve entirely at a higher fluorescence.")
    print("This shows that the entire cell population shifted to a higher fluorescence.")
    print("Conclusion: Antibody B stained ALL cells.\n")

    # Analysis of Histogram C
    print("--- Analysis of Histogram C ---")
    print("The red curve is shifted to the right, but its left tail overlaps significantly with the black curve.")
    print("This indicates the presence of an unstained or weakly stained sub-population.")
    print("Conclusion: Antibody C did NOT stain all cells.\n")

    # Final Conclusion
    print("--- Final Conclusion ---")
    print("Only Histogram B shows a complete shift of the entire population to a higher fluorescence with no remaining unstained cells.")
    print("Therefore, antibody B is the one that stained all cells in the sample.")

analyze_flow_cytometry()