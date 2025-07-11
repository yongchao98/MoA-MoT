def solve_mineral_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Define specimens with their identification and silicate structure classification.
    # A complexity score is assigned: 1 (least complex) to 4 (most complex).
    specimens = [
        {'label': 'A', 'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4}
    ]

    # Sort specimens based on their complexity score.
    sorted_by_complexity = sorted(specimens, key=lambda x: x['complexity'])
    
    # Sort specimens alphabetically for the identification list.
    sorted_by_label = sorted(specimens, key=lambda x: x['label'])

    print("--- Specimen Identification ---")
    for s in sorted_by_label:
        print(f"Specimen {s['label']}: {s['mineral']} ({s['structure']})")
    
    print("\n--- Order of Increasing Silicate Structure Complexity ---")
    
    # Extract labels and structures for the final output string.
    ordered_labels = [s['label'] for s in sorted_by_complexity]
    ordered_structures = [s['structure'] for s in sorted_by_complexity]
    
    # Print the final ordered sequence.
    print(" -> ".join(ordered_labels))
    print(f"({ordered_structures[0]} -> {ordered_structures[1]} -> {ordered_structures[2]} -> {ordered_structures[3]})")


solve_mineral_puzzle()
<<<C, A, B, D>>>