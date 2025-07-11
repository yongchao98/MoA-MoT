def solve_mineral_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1: Identify specimens and classify their silicate structure.
    # The complexity is ranked from 1 (least complex) to 4 (most complex).
    specimens = [
        {'id': 'A', 'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'id': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'id': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'id': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4}
    ]

    # Print the identification of each specimen
    print("Specimen Identification and Classification:")
    for spec in sorted(specimens, key=lambda x: x['id']):
        print(f"Specimen {spec['id']}: {spec['mineral']} ({spec['structure']})")

    # Step 2: Sort the specimens by their structural complexity
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 3: Print the final ordered list
    print("\nOrder of increasing silicate structure complexity:")
    ordered_ids = [spec['id'] for spec in sorted_specimens]
    print(" -> ".join(ordered_ids))
    
    # Store the final answer in the required format
    final_answer = ",".join(ordered_ids)
    print(f"\n<<<C,A,B,D>>>")

solve_mineral_puzzle()