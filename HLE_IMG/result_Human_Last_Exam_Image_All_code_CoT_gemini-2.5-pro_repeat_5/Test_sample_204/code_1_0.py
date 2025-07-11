def solve_mineral_complexity():
    """
    Identifies mineral specimens, classifies their silicate structures,
    and orders them by increasing structural complexity.
    """

    # Step 1: Define the specimens with their identified mineral and silicate structure class.
    specimens = {
        'A': {'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)'},
        'B': {'mineral': 'Mica/Chlorite', 'structure': 'Phyllosilicate (Sheet)'},
        'C': {'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)'},
        'D': {'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)'}
    }

    # Step 2: Define the order of silicate structures from simplest to most complex.
    complexity_order = [
        'Nesosilicate (Isolated Tetrahedra)',
        'Inosilicate (Single Chain)',
        'Phyllosilicate (Sheet)',
        'Tectosilicate (Framework)'
    ]

    # Step 3: Sort the specimens based on the complexity of their silicate structure.
    # The key for sorting is the index of the structure type in the complexity_order list.
    sorted_specimen_keys = sorted(
        specimens.keys(),
        key=lambda k: complexity_order.index(specimens[k]['structure'])
    )

    # Step 4: Print the identification and the final ordered result.
    print("Identifying and ordering specimens by increasing silicate structure complexity:\n")
    for key in sorted_specimen_keys:
        mineral_info = specimens[key]
        print(f"Specimen {key}: {mineral_info['mineral']} is a {mineral_info['structure']}.")

    final_order_str = ', '.join(sorted_specimen_keys)
    print(f"\nThe final order from simplest to most complex structure is: {final_order_str}")

solve_mineral_complexity()
<<<C, A, B, D>>>