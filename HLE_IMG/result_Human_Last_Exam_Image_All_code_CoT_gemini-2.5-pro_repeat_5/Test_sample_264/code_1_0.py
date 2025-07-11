def analyze_flow_cytometry_histograms():
    """
    Analyzes fictional flow cytometry data to determine which antibody stained all cells.
    """

    # Data representation based on visual inspection of the histograms.
    # 'complete' means all cells are stained.
    # 'partial' means only a subpopulation of cells is stained.
    experiments = {
        'A': {
            'staining_pattern': 'partial',
            'reason': 'The red curve shows two peaks. One peak overlaps with the unstained control (black curve), indicating a population of unstained cells remains.'
        },
        'B': {
            'staining_pattern': 'complete',
            'reason': 'The red curve is a single peak, completely shifted to a higher fluorescence compared to the unstained control. This indicates all cells were stained.'
        },
        'C': {
            'staining_pattern': 'partial',
            'reason': 'The red curve shows a major stained peak and a shoulder/small peak that overlaps with the unstained control, indicating some cells remained unstained.'
        }
    }

    print("--- Analysis of Flow Cytometry Histograms ---")
    all_cells_stained_by = []

    for antibody, data in experiments.items():
        print(f"\nAnalysis for Antibody {antibody}:")
        print(f"Result: Staining is {data['staining_pattern']}.")
        print(f"Reasoning: {data['reason']}")
        if data['staining_pattern'] == 'complete':
            all_cells_stained_by.append(antibody)

    print("\n--- Conclusion ---")
    if len(all_cells_stained_by) == 1:
        print(f"Based on the analysis, only antibody '{all_cells_stained_by[0]}' stained all cells in the sample.")
    elif len(all_cells_stained_by) > 1:
        print(f"Based on the analysis, antibodies {', '.join(all_cells_stained_by)} stained all cells in their respective samples.")
    else:
        print("Based on the analysis, none of the antibodies stained all cells in the sample.")

analyze_flow_cytometry_histograms()