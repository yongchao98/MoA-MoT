def solve_mineral_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1: Define the specimens with their identification, silicate class, and complexity level.
    # Complexity is ranked from 1 (least complex) to 4 (most complex).
    specimens = [
        {'id': 'A', 'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'id': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'id': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'id': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4}
    ]

    # Step 2: Sort the specimens based on their complexity level.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 3: Print the results in order.
    print("Ordering specimens by increasing silicate structure complexity:")
    
    order_list = []
    for spec in sorted_specimens:
        print(f"Specimen {spec['id']}: {spec['mineral']}")
        print(f"  - Silicate Structure: {spec['structure']}")
        print(f"  - Complexity arises from how SiO4 tetrahedra are linked. {spec['structure']}s are more complex than the previous and less complex than the next.")
        order_list.append(spec['id'])
        
    final_order = ", ".join(order_list)
    print(f"\nThe final order of increasing complexity is: {final_order}")

solve_mineral_puzzle()