def solve_mineral_puzzle():
    """
    Identifies mineral specimens, classifies their silicate structures,
    and orders them by increasing structural complexity.
    """
    # Step 1-4: Identify specimens and classify their silicate structure.
    # The complexity rank is based on the degree of polymerization of silica tetrahedra:
    # 1: Nesosilicate (isolated)
    # 2: Cyclosilicate (rings)
    # 3: Phyllosilicate (sheets)
    # 4: Tectosilicate (framework)
    specimens = [
        {'label': 'A', 'mineral': 'Tourmaline', 'structure': 'Cyclosilicate', 'complexity_rank': 2},
        {'label': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet Silicate)', 'complexity_rank': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity_rank': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework Silicate)', 'complexity_rank': 4},
    ]

    # Step 5: Sort the specimens by their complexity rank.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity_rank'])

    # Print the identification of each specimen
    print("Specimen Identification and Silicate Structure:")
    for specimen in specimens:
        print(f"  - Specimen {specimen['label']}: {specimen['mineral']} ({specimen['structure']})")
    
    print("\nOrder of increasing silicate structure complexity:")
    
    # Extract the labels in the new sorted order.
    ordered_labels = [s['label'] for s in sorted_specimens]
    
    # Print the final ordered sequence.
    print(" -> ".join(ordered_labels))

solve_mineral_puzzle()