def solve_mineral_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Data for each specimen: Label, Mineral Name, Silicate Structure, Complexity Rank
    specimens = [
        {'label': 'A', 'mineral': 'Pyroxene', 'structure': 'Single-chain silicate (Inosilicate)', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica', 'structure': 'Sheet silicate (Phyllosilicate)', 'complexity': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Isolated tetrahedra (Nesosilicate)', 'complexity': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Framework silicate (Tectosilicate)', 'complexity': 4}
    ]

    print("--- Specimen Identification ---")
    # Sort by label for clear presentation
    specimens_by_label = sorted(specimens, key=lambda x: x['label'])
    for s in specimens_by_label:
        print(f"Specimen {s['label']}: {s['mineral']} ({s['structure']})")

    # Sort specimens by complexity rank
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    print("\n--- Order of Increasing Silicate Complexity ---")
    ordered_labels = [s['label'] for s in sorted_specimens]
    
    # Print the final equation/order
    print("The order from least complex to most complex is:")
    print(" -> ".join(ordered_labels))

    # Final answer in the required format
    final_answer = ", ".join(ordered_labels)
    print(f"\n<<< {final_answer} >>>")

solve_mineral_puzzle()