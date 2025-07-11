def solve_mineral_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    specimens = [
        {'label': 'A', 'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4},
    ]

    print("Specimen Identification:")
    # Print identification for each specimen based on their labels A, B, C, D
    for label_char in sorted([s['label'] for s in specimens]):
        for specimen in specimens:
            if specimen['label'] == label_char:
                print(f"- Specimen {specimen['label']}: {specimen['mineral']} ({specimen['structure']})")
                break
    
    # Sort specimens by complexity
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])
    
    # Get the ordered labels
    ordered_labels = [s['label'] for s in sorted_specimens]
    
    print("\nOrder of specimens by increasing silicate structure complexity:")
    # The final equation is the sequence of labels
    final_equation = " < ".join(ordered_labels)
    print(final_equation)

solve_mineral_puzzle()
