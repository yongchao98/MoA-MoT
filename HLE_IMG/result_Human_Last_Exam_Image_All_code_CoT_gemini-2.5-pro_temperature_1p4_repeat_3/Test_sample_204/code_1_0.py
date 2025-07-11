def solve_silicate_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    specimens = [
        {'label': 'A', 'mineral': 'Pyroxene', 'structure': 'Single-Chain Silicate (Inosilicate)', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica/Chlorite', 'structure': 'Sheet Silicate (Phyllosilicate)', 'complexity': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Isolated Tetrahedra (Nesosilicate)', 'complexity': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Framework Silicate (Tectosilicate)', 'complexity': 4}
    ]

    print("Specimen Identification:")
    # Sort alphabetically to present initial identification clearly
    for specimen in sorted(specimens, key=lambda x: x['label']):
        print(f"Specimen {specimen['label']}: {specimen['mineral']} ({specimen['structure']})")

    # Sort by complexity
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    print("\nOrder of increasing silicate structure complexity:")
    
    # Get the labels in the correct order
    ordered_labels = [s['label'] for s in sorted_specimens]
    
    # Print the final ordered list
    print(" -> ".join(ordered_labels))
    
    # Prepare the final answer string
    final_answer = ", ".join(ordered_labels)
    print(f"\n<<<C, A, B, D>>>")

solve_silicate_puzzle()