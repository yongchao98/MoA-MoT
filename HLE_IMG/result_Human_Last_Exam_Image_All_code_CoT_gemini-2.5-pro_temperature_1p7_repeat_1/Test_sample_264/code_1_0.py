def analyze_flow_cytometry_histograms():
    """
    This function analyzes the provided flow cytometry data based on visual interpretation
    to determine which antibodies stained all cells in their respective samples.
    """

    print("--- Analysis of Flow Cytometry Histograms ---")
    print("\nPrinciple:")
    print("If an antibody stains ALL cells, the entire population's fluorescence increases.")
    print("On the histogram, this means the red curve (stained) completely shifts to the right,")
    print("leaving no cell population behind at the original, low fluorescence of the black curve (unstained).\n")

    # Dictionary to store the analysis result for each antibody
    # True means it stained all cells, False means it stained a subset.
    analysis = {}

    # Step 1: Analyze Histogram A
    print("--- Analyzing Antibody A ---")
    print("Observation: The red curve is a single peak, shifted to the right of the black control curve. No peak remains at the original 'unstained' position.")
    print("Conclusion: Antibody A stained ALL cells.")
    analysis['A'] = True

    # Step 2: Analyze Histogram B
    print("\n--- Analyzing Antibody B ---")
    print("Observation: The red curve shows two peaks. One peak is on the far right (stained cells), but another peak overlaps perfectly with the black control curve (unstained cells).")
    print("Conclusion: Antibody B stained only a SUBSET of cells.")
    analysis['B'] = False

    # Step 3: Analyze Histogram C
    print("\n--- Analyzing Antibody C ---")
    print("Observation: The red curve is a single peak, completely shifted to the right of the black control curve, with a clear separation.")
    print("Conclusion: Antibody C stained ALL cells.")
    analysis['C'] = True

    # Final Summary
    print("\n--- Final Summary ---")
    all_cells_stained_by = [antibody for antibody, result in analysis.items() if result]
    
    if all_cells_stained_by:
        result_string = " and ".join(all_cells_stained_by)
        print(f"The antibodies that stained all cells are: {result_string}.")
        # Matching with the answer choices
        if sorted(all_cells_stained_by) == ['A', 'C']:
            answer_choice = "C"
            print(f"This corresponds to answer choice {answer_choice}.")
        else:
            print("The result does not match the provided answer choices.")
    else:
        print("No antibody stained all the cells.")


# Run the analysis
analyze_flow_cytometry_histograms()