def solve_mineral_ordering():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Define the specimens with their identified mineral and silicate structure
    specimens = [
        {'id': 'A', 'mineral': 'Pyroxene (e.g., Augite)', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'id': 'B', 'mineral': 'Mica/Chlorite', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'id': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'id': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4}
    ]

    # Sort specimens based on the complexity rank
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Print the identification of each specimen
    print("Specimen Identification:")
    for spec in sorted(specimens, key=lambda x: x['id']):
        print(f"Specimen {spec['id']}: {spec['mineral']} ({spec['structure']})")

    # Print the final ordered list
    print("\nSpecimens ordered by increasing silicate structure complexity:")
    order_string = " -> ".join([spec['id'] for spec in sorted_specimens])
    print(order_string)

solve_mineral_ordering()
<<<C -> A -> B -> D>>>