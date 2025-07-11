def solve_mineral_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1: Define the specimens with their identified mineral type and silicate structure.
    # We assign a complexity score where a lower number is less complex.
    # Nesosilicate (1) < Cyclosilicate (2) < Phyllosilicate (3) < Tectosilicate (4)
    specimens = [
        {'label': 'A', 'mineral': 'Tourmaline', 'structure': 'Cyclosilicate', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate', 'complexity': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate', 'complexity': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate', 'complexity': 4}
    ]

    # Step 2: Sort the specimens based on their complexity score.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 3: Print the identification for each specimen.
    print("Specimen Identification:")
    # Create a copy sorted by label for clear printing
    specimens_by_label = sorted(specimens, key=lambda x: x['label'])
    for spec in specimens_by_label:
        print(f"Specimen {spec['label']}: {spec['mineral']} ({spec['structure']})")

    # Step 4: Print the final ordered list.
    print("\nOrder of increasing silicate structure complexity:")
    
    # Extract the labels from the sorted list
    ordered_labels = [spec['label'] for spec in sorted_specimens]
    
    # Build the final equation string showing the order
    final_order_string = " -> ".join(ordered_labels)
    
    print(f"{sorted_specimens[0]['label']} ({sorted_specimens[0]['structure']}) -> {sorted_specimens[1]['label']} ({sorted_specimens[1]['structure']}) -> {sorted_specimens[2]['label']} ({sorted_specimens[2]['structure']}) -> {sorted_specimens[3]['label']} ({sorted_specimens[3]['structure']})")
    print(f"\nThe final order is: {final_order_string}")

solve_mineral_puzzle()