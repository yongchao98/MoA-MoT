def analyze_flow_cytometry():
    """
    Analyzes conceptual data from flow cytometry histograms to determine which antibody stained all cells.
    """
    # Step 1: Define the interpretation of each histogram based on visual analysis.
    # 'partial': Indicates that only a subset of cells was stained, resulting in two populations
    #            (stained and unstained) in the antibody-treated sample.
    # 'all': Indicates that all cells were stained, resulting in the entire population shifting
    #        to a higher fluorescence.
    staining_results = {
        'A': 'partial',
        'B': 'all',
        'C': 'partial'
    }

    # Step 2: Provide a detailed explanation of the principles.
    print("--- Understanding the Histograms ---")
    print("In flow cytometry, we measure the fluorescence of individual cells.")
    print(" - The black curve shows the 'unstained' control population, which has low baseline fluorescence.")
    print(" - The red curve shows the 'antibody-stained' population.")
    print("\nIf an antibody stains ALL cells:")
    print("The entire red curve will shift to a higher fluorescence (to the right) compared to the black curve.")
    print("There will be no significant population of cells left at the low, unstained fluorescence level.")
    print("\nIf an antibody stains only a SUBSET of cells:")
    print("The red curve will show two populations: one that remains at the low, unstained fluorescence level and one that shifts to a higher fluorescence.")

    # Step 3: Analyze each case and find the solution.
    print("\n--- Analysis of Each Antibody ---")
    all_cells_stained_by = []
    for antibody, result in staining_results.items():
        if result == 'all':
            print(f"Antibody {antibody}: The red curve is completely shifted to the right. This indicates that ALL cells were stained.")
            all_cells_stained_by.append(antibody)
        else:
            print(f"Antibody {antibody}: The red curve shows an overlap with the black curve. This indicates that only a SUBSET of cells was stained.")

    # Step 4: Print the final conclusion.
    print("\n--- Conclusion ---")
    if len(all_cells_stained_by) == 1:
        print(f"The only antibody that stained all cells in the sample is Antibody {all_cells_stained_by[0]}.")
    else:
        print(f"The antibodies that stained all cells are: {', '.join(all_cells_stained_by)}.")

analyze_flow_cytometry()
<<<E>>>