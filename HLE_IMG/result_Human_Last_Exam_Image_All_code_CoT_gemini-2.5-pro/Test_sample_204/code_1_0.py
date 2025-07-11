def solve_mineral_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1 & 2: Define each specimen with its identification, silicate group, and a complexity score.
    # The complexity score is based on the degree of polymerization of silica tetrahedra:
    # 1: Nesosilicate (isolated)
    # 2: Inosilicate (chain)
    # 3: Phyllosilicate (sheet)
    # 4: Tectosilicate (framework)
    specimens = [
        {'label': 'C', 'mineral': 'Olivine', 'group': 'Nesosilicate', 'complexity': 1},
        {'label': 'A', 'mineral': 'Pyroxene', 'group': 'Inosilicate (single chain)', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica', 'group': 'Phyllosilicate (sheet)', 'complexity': 3},
        {'label': 'D', 'mineral': 'Quartz', 'group': 'Tectosilicate (framework)', 'complexity': 4}
    ]

    # Step 3: Sort the specimens based on their assigned complexity score.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 4: Print the identification of each specimen and the final ordered list.
    print("Specimen Identification and Silicate Group:")
    # Create a correctly ordered dictionary for printing identifications
    id_dict = {item['label']: item for item in specimens}
    for label in sorted(id_dict.keys()):
        item = id_dict[label]
        print(f"Specimen {item['label']}: {item['mineral']} ({item['group']})")
    
    print("\nOrder of increasing silicate structure complexity:")
    
    # Build the descriptive output string
    order_description = []
    for specimen in sorted_specimens:
        order_description.append(f"{specimen['label']} ({specimen['mineral']})")
    
    print(" -> ".join(order_description))

    # Extract the final ordered labels
    final_order_labels = [specimen['label'] for specimen in sorted_specimens]
    
    print("\nFinal Answer (the ordered sequence):")
    print(", ".join(final_order_labels))

solve_mineral_puzzle()